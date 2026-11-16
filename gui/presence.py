import time

CLIENT_ID = "1571825555970588744"


class Presence:
    """Rich Presence Discord. Silencieux si Discord est indisponible."""

    def __init__(self, client_id=CLIENT_ID):
        self.rpc = None
        self.debut = int(time.time())
        try:
            from pypresence import Presence as RPC
            rpc = RPC(client_id)
            rpc.connect()
            self.rpc = rpc
            self.maj("En cours de lancement")
        except Exception:
            self.rpc = None

    def maj(self, etat):
        if self.rpc is None:
            return
        try:
            self.rpc.update(details="Bio-informatique", state=etat,
                            start=self.debut,
                            large_image="logo", large_text="bioapp")
        except Exception:
            self.rpc = None

    def fermer(self):
        if self.rpc is not None:
            try:
                self.rpc.close()
            except Exception:
                pass
