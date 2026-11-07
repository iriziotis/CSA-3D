.PHONY: help install install-dev test lint format type-check clean docs

help:
	@echo "CSA-3D Development Tools"
	@echo "======================="
	@echo ""
	@echo "Available targets:"
	@echo "  install          Install package in production mode"
	@echo "  install-dev      Install with development dependencies"
	@echo "  test             Run tests with pytest"
	@echo "  test-cov         Run tests with coverage report"
	@echo "  lint             Run flake8, isort, and black checks"
	@echo "  format           Auto-format code with black and isort"
	@echo "  type-check       Run mypy type checking"
	@echo "  clean            Remove build artifacts and cache files"
	@echo "  docs             Build HTML documentation"
	@echo ""

install:
	pip install -e .

install-dev:
	pip install -e ".[dev,docs,examples]"

test:
	pytest tests/ -v

test-cov:
	pytest tests/ -v --cov=csa3d --cov-report=html --cov-report=term-missing

lint:
	@echo "Running isort check..."
	isort --check-only csa3d tests
	@echo "Running black check..."
	black --check csa3d tests
	@echo "Running flake8..."
	flake8 csa3d tests
	@echo "All checks passed!"

format:
	@echo "Running isort..."
	isort csa3d tests
	@echo "Running black..."
	black csa3d tests
	@echo "Code formatted!"

type-check:
	mypy csa3d --ignore-missing-imports

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .mypy_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .coverage -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name htmlcov -exec rm -rf {} + 2>/dev/null || true
	rm -rf build/ dist/ .eggs/ *.egg-info/

docs:
	@echo "Building documentation..."
	cd docs && make html
	@echo "Documentation built in docs/_build/html/"
