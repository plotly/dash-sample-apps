from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="pkg-NSKHALDI", version="1.0", author="Anass Khaldi", packages=find_packages()
)
