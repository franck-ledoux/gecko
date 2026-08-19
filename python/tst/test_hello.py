"""Smoke test for the gecko Python extension — see docs/user-guide/python.md."""

import gecko


def test_hello():
    assert gecko.hello() == "Hello from gecko!"
