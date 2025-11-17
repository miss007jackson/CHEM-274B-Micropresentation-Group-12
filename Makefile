# Makefile for ProteinGraph project with virtual environment

PYTHON=python3
VENV=venv
SCRIPT=protein.py
PDB_FILE=1ubq.pdb
PIP=$(VENV)/bin/pip
PYTHON_VENV=$(VENV)/bin/python

.PHONY: all venv install run clean

# Default target
all: venv install run

# Create virtual environment
venv:
	$(PYTHON) -m venv $(VENV)

# Check Graphviz system dependencies (for pygraphviz)
check-graphviz:
	@if ! dpkg -s graphviz >/dev/null 2>&1; then \
		echo "ERROR: Graphviz is not installed."; \
		echo "Install it using: sudo apt install graphviz graphviz-dev"; \
		exit 1; \
	fi

# Install required packages
install: venv
	$(PIP) install --upgrade pip
	$(PIP) install biopython networkx matplotlib pygraphviz

# Run the script inside the virtual environment
run: install
	$(PYTHON_VENV) $(SCRIPT)

# Clean generated files and environment
clean:
	rm -rf $(VENV) *.png __pycache__