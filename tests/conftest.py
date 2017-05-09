import os
import random
import string
import pytest
import subprocess


def pytest_report_header(config):
    freeze = subprocess.check_output(['pip', 'freeze'])
    return "QuakeML test runner: qmlutil" + '\n' + freeze


def pytest_addoption(parser):
    parser.addoption("--integration", action="store_true", 
        help="run integration tests")
    parser.addoption("--writefiles", action="store_true", 
        help="Write test files for inspection")


@pytest.fixture
def random_filename(tmpdir_factory):
    def make_random_filename(ext=''):
        dir = str(tmpdir_factory.mktemp('qmlutil').realpath())
        fname = ''.join(random.choice(string.ascii_lowercase)
                        for _ in range(10))
        return os.path.join(dir, fname + ext)
    return make_random_filename