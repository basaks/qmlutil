from os.path import isfile
import pytest
from obspy.core.event import read_events
try:
    from obspy.io.seiscomp.event import _write_sc3ml
    SKIP = False
except:
    SKIP = True


@pytest.mark.skipif(SKIP, reason='Installed version of obspy does not support'
                                 'conversion to seiscompML')
def test_functionality(random_filename):
    catalog = read_events()
    sc3xml = random_filename(ext='.xml')
    _write_sc3ml(catalog=catalog, filename=sc3xml)
    assert isfile(sc3xml)




