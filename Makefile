.PHONY: verify test

verify:
	python3 scripts/verify_workshop.py

test:
	python3 -m unittest discover -s tests -v
