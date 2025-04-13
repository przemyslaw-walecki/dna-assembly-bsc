import time

class Authenticator:
    def __init__(self, user_db):
        self.user_db = user_db

    def login(self, username, password):
        start = time.time()
        if username in self.user_db and self.user_db[username] == password:
            return True, time.time() - start
        return False, time.time() - start
