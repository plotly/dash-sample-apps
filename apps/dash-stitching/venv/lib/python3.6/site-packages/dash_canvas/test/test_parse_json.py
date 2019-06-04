from dash_canvas.utils import parse_jsonfile


def test_parse_jsonfile():
    shape = (433, 640)
    mask = parse_jsonfile('tests/data_test.json', shape=shape)
    assert mask.sum() > 0
    assert mask.shape == shape
