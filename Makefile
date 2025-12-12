# Prevent running inside conda
ifneq ($(CONDA_DEFAULT_ENV),)
$(error Do NOT run make inside a conda environment)
endif

# Use preferred Python
PYTHON := python3
VENV_DIR := .venv
PIP := $(VENV_DIR)/bin/pip
PY := $(VENV_DIR)/bin/python

# Default: set up env
.PHONY: all
all: venv install

# Try to create virtual environment; if it fails, continue without it
.PHONY: venv
venv:
	@if [ ! -d "$(VENV_DIR)" ]; then \
		echo ">>> Attempting to create virtual environment in $(VENV_DIR)"; \
		if $(PYTHON) -m venv $(VENV_DIR) 2>/dev/null; then \
			echo ">>> Virtual environment created successfully"; \
		else \
			echo ">>> WARNING: Failed to create virtual environment. Falling back to system Python."; \
		fi \
	else \
		echo ">>> Virtual environment already exists"; \
	fi

# Install dependencies (with user confirmation)
.PHONY: install
install: venv
	@echo ">>> Do you want to install dependencies from requirements.txt? (y/n)"
	@read ans; \
	if [ "$$ans" != "y" ]; then \
		echo ">>> Installation cancelled by user"; \
		exit 1; \
	fi; \
	echo ">>> Installing packages..."; \
	if [ -x "$(PIP)" ]; then \
		echo ">>> Using virtual environment pip"; \
		$(PIP) install --upgrade pip; \
		$(PIP) install -r requirements.txt; \
	else \
		echo ">>> Using system pip (no venv available)"; \
		$(PYTHON) -m pip install --user --upgrade pip; \
		$(PYTHON) -m pip install --user -r requirements.txt; \
	fi


# Run main script
.PHONY: run
run: install
	@echo ">>> Running main program"
	@if [ -x "$(PY)" ]; then \
		echo ">>> Using virtual environment Python"; \
		$(PY) main.py; \
	else \
		echo ">>> Using system Python"; \
		$(PYTHON) main.py; \
	fi

# Remove the virtual env
.PHONY: clean
clean:
	@echo ">>> Removing virtual environment"
	rm -rf $(VENV_DIR)
