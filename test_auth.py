import unittest
from auth import Authenticator


class Authenticator:
    pass

class TestAuthenticator(unittest.TestCase):
    def setUp(self):
        self.auth = Authenticator({"user": "pass"})

    def test_login_success(self):
        success, _ = self.auth.login("user", "pass")
        self.assertTrue(success)

    def test_login_performance(self):
        _, duration = self.auth.login("user", "pass")
        self.assertLess(duration, 0.5)

if __name__ == '__main__':
    unittest.main()
