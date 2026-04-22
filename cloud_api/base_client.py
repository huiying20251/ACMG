#!/usr/bin/env python3
"""
Base HTTP client for cloud API calls.
Provides common retry logic, timeout handling, and error processing.
"""

import time
import logging
from typing import Optional, Dict, Any
from dataclasses import dataclass

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)


@dataclass
class APIResponse:
    """Standardized API response wrapper."""
    success: bool
    data: Optional[Any] = None
    error: Optional[str] = None
    status_code: Optional[int] = None
    source: str = "api"  # "api" or "cache" or "local"


class BaseCloudClient:
    """
    Base class for cloud API clients.
    Provides common HTTP functionality with retry logic and timeout handling.
    """

    DEFAULT_TIMEOUT = 30  # seconds
    DEFAULT_RETRIES = 3
    DEFAULT_BACKOFF_FACTOR = 0.5

    def __init__(
        self,
        base_url: str,
        timeout: int = DEFAULT_TIMEOUT,
        max_retries: int = DEFAULT_RETRIES,
        backoff_factor: float = DEFAULT_BACKOFF_FACTOR,
        headers: Optional[Dict[str, str]] = None,
    ):
        """
        Initialize base cloud client.

        Args:
            base_url: Base URL for API endpoint
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
            backoff_factor: Backoff factor for retries
            headers: Default headers to include in requests
        """
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.session = self._create_session(max_retries, backoff_factor)
        self.headers = headers or {}

    def _create_session(self, max_retries: int, backoff_factor: float) -> requests.Session:
        """
        Create a requests session with retry logic.

        Args:
            max_retries: Maximum number of retry attempts
            backoff_factor: Backoff factor for retries

        Returns:
            Configured requests session
        """
        session = requests.Session()

        retry_strategy = Retry(
            total=max_retries,
            backoff_factor=backoff_factor,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET", "POST"],
        )

        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("http://", adapter)
        session.mount("https://", adapter)

        return session

    def _make_request(
        self,
        method: str,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: Optional[int] = None,
    ) -> APIResponse:
        """
        Make an HTTP request with error handling.

        Args:
            method: HTTP method (GET, POST, etc.)
            endpoint: API endpoint path
            params: Query parameters
            data: Request body data
            headers: Additional headers
            timeout: Request timeout override

        Returns:
            APIResponse with success status and data/error
        """
        url = f"{self.base_url}/{endpoint.lstrip('/')}"
        request_headers = {**self.headers, **(headers or {})}
        request_timeout = timeout or self.timeout

        try:
            response = self.session.request(
                method=method,
                url=url,
                params=params,
                json=data,
                headers=request_headers,
                timeout=request_timeout,
            )

            response.raise_for_status()

            # Try to parse JSON
            try:
                response_data = response.json()
            except ValueError:
                response_data = response.text

            return APIResponse(
                success=True,
                data=response_data,
                status_code=response.status_code,
            )

        except requests.exceptions.Timeout:
            logger.error(f"Request timeout for {url}")
            return APIResponse(
                success=False,
                error=f"Request timeout after {request_timeout}s",
                status_code=408,
            )

        except requests.exceptions.ConnectionError as e:
            logger.error(f"Connection error for {url}: {e}")
            return APIResponse(
                success=False,
                error=f"Connection error: {str(e)}",
                status_code=503,
            )

        except requests.exceptions.HTTPError as e:
            logger.error(f"HTTP error for {url}: {e}")
            return APIResponse(
                success=False,
                error=str(e),
                status_code=response.status_code,
            )

        except Exception as e:
            logger.error(f"Unexpected error for {url}: {e}")
            return APIResponse(
                success=False,
                error=f"Unexpected error: {str(e)}",
                status_code=500,
            )

    def get(
        self,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: Optional[int] = None,
    ) -> APIResponse:
        """
        Make a GET request.

        Args:
            endpoint: API endpoint path
            params: Query parameters
            headers: Additional headers
            timeout: Request timeout override

        Returns:
            APIResponse with success status and data/error
        """
        return self._make_request("GET", endpoint, params=params, headers=headers, timeout=timeout)

    def post(
        self,
        endpoint: str,
        data: Optional[Dict[str, Any]] = None,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: Optional[int] = None,
    ) -> APIResponse:
        """
        Make a POST request.

        Args:
            endpoint: API endpoint path
            data: Request body data
            params: Query parameters
            headers: Additional headers
            timeout: Request timeout override

        Returns:
            APIResponse with success status and data/error
        """
        return self._make_request("POST", endpoint, params=params, data=data, headers=headers, timeout=timeout)


class RateLimiter:
    """
    Simple rate limiter for API calls.
    """

    def __init__(self, calls_per_second: float = 10.0):
        """
        Initialize rate limiter.

        Args:
            calls_per_second: Maximum API calls per second
        """
        self.calls_per_second = calls_per_second
        self.min_interval = 1.0 / calls_per_second
        self.last_call_time = 0.0

    def wait(self) -> None:
        """Wait if necessary to respect rate limits."""
        elapsed = time.time() - self.last_call_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_call_time = time.time()
