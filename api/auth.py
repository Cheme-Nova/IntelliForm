import os
from typing import Optional

import requests


SUPABASE_URL = os.getenv("SUPABASE_URL", "").rstrip("/")
SUPABASE_ANON_KEY = os.getenv("SUPABASE_ANON_KEY", "")
REQUIRE_SIGNIN = os.getenv("INTELLIFORM_REQUIRE_SIGNIN", "1").lower() not in {"0", "false", "no"}


def is_auth_enabled() -> bool:
    return bool(SUPABASE_URL and SUPABASE_ANON_KEY)


def extract_bearer_token(request) -> str:
    header = request.headers.get("authorization", "")
    if header.lower().startswith("bearer "):
        return header[7:].strip()
    return ""


def verify_supabase_user(access_token: str) -> Optional[dict]:
    if not access_token or not is_auth_enabled():
        return None
    try:
        response = requests.get(
            f"{SUPABASE_URL}/auth/v1/user",
            headers={
                "apikey": SUPABASE_ANON_KEY,
                "Authorization": f"Bearer {access_token}",
            },
            timeout=10,
        )
        if response.status_code != 200:
            return None
        return response.json()
    except Exception:
        return None
