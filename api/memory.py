from datetime import datetime

class MemoryBus:
    def __init__(self):
        self.store = []

    def record(self, event_type, payload, vertical=None):
        self.store.append({
            "type": event_type,
            "payload": payload,
            "vertical": vertical,
            "timestamp": datetime.utcnow().isoformat()
        })

    def query(self, filter_fn):
        return [m for m in self.store if filter_fn(m)]

    def recent(self, n=10):
        return self.store[-n:]

memory = MemoryBus()
