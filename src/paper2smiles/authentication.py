"""Verify Portal gateway identities before any job access."""

import base64
import hashlib
import hmac
import json
import time

from starlette.exceptions import HTTPException


class GatewayAuthenticator:
    """Verify Portal gateway identities before any job access.

    Parameters
    ----------
    config_provider : callable
        Return the signing secret and expected gateway audience.
    """

    def __init__(self, config_provider):
        """Bind the service configuration provider."""
        self.settings = config_provider

    def owner(self, request):
        """Verify the gateway token and recover its authenticated subject.

        Parameters
        ----------
        request : starlette.requests.Request
            Incoming authenticated Portal gateway request.

        Returns
        -------
        str
            Owner identifier bound to this request method and path.

        Raises
        ------
        starlette.exceptions.HTTPException
            If the signature, subject, expiry, method, or path is invalid.
        """
        token = request.headers.get("authorization", "").removeprefix(
            "Bearer "
        )
        secret = self.settings()["identity_secret"]
        try:
            header, payload, signature = token.split(".")

            def decode(s):
                return base64.urlsafe_b64decode(s + "=" * (-len(s) % 4))

            if len(secret) < 32:
                raise ValueError()
            if json.loads(decode(header)).get("alg") != "HS256":
                raise ValueError()
            expected = hmac.new(
                secret.encode(),
                (header + "." + payload).encode(),
                hashlib.sha256,
            ).digest()
            if not hmac.compare_digest(expected, decode(signature)):
                raise ValueError()
            claims = json.loads(decode(payload))
            now = time.time()
            if (
                claims.get("iss") != "aidd-portal"
                or claims.get("aud") != self.settings()["audience"]
                or claims.get("exp", 0) <= now
                or (claims.get("iat", now + 100) > now + 10)
            ):
                raise ValueError()
            if (
                claims.get("path") != request.url.path
                or claims.get("method") != request.method
                or (not isinstance(claims.get("sub"), str))
                or (not claims["sub"])
            ):
                raise ValueError()
            return claims["sub"]
        except (ValueError, TypeError, KeyError, json.JSONDecodeError):
            raise HTTPException(401, "Authentication required")
