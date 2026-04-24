import os
import time
from collections import deque
from typing import Optional, Tuple


FREE_TIER_ENABLED = os.getenv("INTELLIFORM_FREE_TIER", "1").lower() not in {"0", "false", "no"}
FREE_TIER_MAX_REQUESTS_PER_HOUR = int(os.getenv("INTELLIFORM_FREE_TIER_MAX_REQUESTS_PER_HOUR", "8"))
FREE_TIER_MAX_QSAR_PER_HOUR = int(os.getenv("INTELLIFORM_FREE_TIER_MAX_QSAR_PER_HOUR", "20"))
FREE_TIER_REQUIRE_API_KEY = os.getenv("INTELLIFORM_FREE_TIER_REQUIRE_API_KEY", "0").lower() in {"1", "true", "yes"}
PUBLIC_API_KEY = os.getenv("INTELLIFORM_PUBLIC_API_KEY", "")


_request_windows = {}


def _window_key(client_id: str, bucket: str) -> tuple[str, str]:
    return client_id, bucket


def _allow(client_id: str, bucket: str, limit: int, window_seconds: int = 3600) -> tuple[bool, int]:
    now = time.time()
    key = _window_key(client_id, bucket)
    queue = _request_windows.setdefault(key, deque())

    while queue and (now - queue[0]) >= window_seconds:
        queue.popleft()

    if len(queue) >= limit:
        retry_after = max(1, int(window_seconds - (now - queue[0])))
        return False, retry_after

    queue.append(now)
    return True, 0


def get_client_id(request) -> str:
    forwarded = request.headers.get("x-forwarded-for", "")
    if forwarded:
        return forwarded.split(",")[0].strip()
    if request.client and request.client.host:
        return request.client.host
    return "unknown"


def validate_public_access(request, bucket: str = "formulate") -> Tuple[bool, int, Optional[str]]:
    if not FREE_TIER_ENABLED:
        return True, 0, None

    if FREE_TIER_REQUIRE_API_KEY:
        provided = (request.headers.get("x-api-key") or "").strip()
        if not PUBLIC_API_KEY or provided != PUBLIC_API_KEY:
            return False, 0, "Missing or invalid public API key."

    client_id = get_client_id(request)
    limit = FREE_TIER_MAX_QSAR_PER_HOUR if bucket == "qsar" else FREE_TIER_MAX_REQUESTS_PER_HOUR
    allowed, retry_after = _allow(client_id, bucket, limit)
    if not allowed:
        return False, retry_after, f"Free-tier limit reached for {bucket}. Try again later."
    return True, 0, None
