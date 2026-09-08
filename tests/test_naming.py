import urllib.error
from unittest.mock import MagicMock, patch

from core.naming import get_iupac_name


def _mock_response(text: str):
    response = MagicMock()
    response.__enter__.return_value.read.return_value = text.encode("utf-8")
    return response


def test_get_iupac_name_returns_name_on_success():
    with patch("core.naming.urllib.request.urlopen", return_value=_mock_response("ethanol\n")):
        assert get_iupac_name("CCO") == "ethanol"


def test_get_iupac_name_returns_none_on_404():
    error = urllib.error.HTTPError(url="", code=404, msg="Not Found", hdrs=None, fp=None)
    with patch("core.naming.urllib.request.urlopen", side_effect=error):
        assert get_iupac_name("some_unknown_structure") is None


def test_get_iupac_name_returns_none_on_other_http_error():
    error = urllib.error.HTTPError(url="", code=400, msg="Bad Request", hdrs=None, fp=None)
    with patch("core.naming.urllib.request.urlopen", side_effect=error):
        assert get_iupac_name("CCO") is None


def test_get_iupac_name_returns_none_on_network_error():
    with patch("core.naming.urllib.request.urlopen", side_effect=urllib.error.URLError("no network")):
        assert get_iupac_name("CCO") is None


def test_get_iupac_name_returns_none_on_unexpected_exception():
    with patch("core.naming.urllib.request.urlopen", side_effect=TimeoutError("timed out")):
        assert get_iupac_name("CCO") is None


def test_get_iupac_name_returns_none_for_empty_response():
    with patch("core.naming.urllib.request.urlopen", return_value=_mock_response("  \n")):
        assert get_iupac_name("CCO") is None


def test_get_iupac_name_does_not_raise_and_never_hits_network_when_mocked():
    # Sanity check that the function never actually opens a socket in tests -
    # if urlopen weren't mocked, this would attempt a real HTTP request.
    with patch("core.naming.urllib.request.urlopen") as mock_urlopen:
        mock_urlopen.return_value = _mock_response("benzene")
        result = get_iupac_name("c1ccccc1")

        assert result == "benzene"
        assert mock_urlopen.called
