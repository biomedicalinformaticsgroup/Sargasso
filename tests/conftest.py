import pytest
from pkg_resources import resource_filename


@pytest.fixture(scope="module")
def test_files():
    file = resource_filename('ddp2neo4j.resources.obo', 'hpo.obo')
