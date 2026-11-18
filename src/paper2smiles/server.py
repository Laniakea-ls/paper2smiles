"""Compose the extraction API from independently testable service objects."""

from collections.abc import Callable

from starlette.applications import Starlette
from starlette.exceptions import HTTPException
from starlette.routing import Route

from paper2smiles.api import ExtractionAPI
from paper2smiles.authentication import GatewayAuthenticator
from paper2smiles.jobs.presentation import JobPresenter
from paper2smiles.jobs.repository import JobRepository
from paper2smiles.jobs.runner import JobRunner
from paper2smiles.settings import settings


def create_app(config_provider: Callable[[], dict] = settings) -> Starlette:
    """Create an API application with its own job-service dependencies.

    Parameters
    ----------
    config_provider : callable, optional
        Return service configuration. Tests may supply an isolated provider.

    Returns
    -------
    starlette.applications.Starlette
        Application with owner-scoped routes and a supervised queue worker.
    """
    repository = JobRepository(config_provider)
    authenticator = GatewayAuthenticator(config_provider)
    presenter = JobPresenter(repository)
    api = ExtractionAPI(repository, authenticator, presenter)
    runner = JobRunner(repository)
    app = Starlette(
        routes=[
            Route("/health", api.health),
            Route("/history", api.history),
            Route("/jobs", api.jobs, methods=["GET", "POST"]),
            Route("/jobs/{job_id}", api.detail),
            Route("/jobs/{job_id}/cancel", api.cancel, methods=["POST"]),
            Route("/jobs/{job_id}/retry", api.retry, methods=["POST"]),
            Route("/jobs/{job_id}/files/{path:path}", api.artifact),
        ],
        lifespan=runner.lifespan,
        exception_handlers={HTTPException: api.error},
    )
    app.state.repository = repository
    app.state.runner = runner
    return app


app = create_app()
