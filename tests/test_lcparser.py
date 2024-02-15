from lcparser.lcparser import LCData, parse_lc_textfile

metadata = {'File Path': 'fake/file/path.txt',
            'Channel': 'channel',
            'Injection Information': {'Injection Volume (µL)': '25.0'},
            'Chromatogram Data Information': {'Units': ['Time (min)', 'Step (s)', 'Value (EU)']},
            'Signal Parameter Information': {}}

data = [[0, '0.5', 0],
        [1, '0.5', 2],
        [2, '0.5', 6],
        [3, '0.5', 1],
        [4, '0.5', 0],
        [5, '0.5', 1],
        [6, '0.5', 1],
        [7, '0.5', 7],
        [8, '0.5', 9],
        [9, '0.5', 3],
        [10, '0.5', 0]]

test_obj = LCData(metadata, data)


def test_parse_lc_textfile():
    lcdata_obj = parse_lc_textfile("data/IgG Vtag 1_ACQUITY FLR ChA.txt")
    assert type(lcdata_obj) == LCData
    assert (type(lcdata_obj.metadata)) == dict
    assert (type(lcdata_obj.data)) == list


def test_getters():
    assert test_obj.get_time() == [0,1,2,3,4,5,6,7,8,9,10]
    assert test_obj.get_step() == ['0.5','0.5','0.5','0.5','0.5','0.5','0.5','0.5','0.5','0.5','0.5']
    assert test_obj.get_value() == [0,2,6,1,0,1,1,7,9,3,0]


def test_process_chromatogram():
    # This throws a warning about np.matrix, but it is from the scipy/baseline stuff.
    assert not hasattr(test_obj, 'ydata_processed')
    test_obj.process_raw_chromatogram()
    assert hasattr(test_obj, 'ydata_processed')


def test_detect_peaks():
    assert not hasattr(test_obj, 'peaks')
    test_obj.detect_peaks()
    assert len(test_obj.peaks.get('peak_index')) == 2
    assert test_obj.peaks.get('peak_index')[0] == 2
    assert test_obj.peaks.get('peak_index')[1] == 8
    assert 'peak_props' in test_obj.peaks.keys()


def test_integrate_peaks():
    test_obj.integrate_peaks()
    assert 'peak_areas' in test_obj.peaks.keys()
    assert len(test_obj.peaks["peak_areas"]) == 2


def test_adjust_peak_bases():
    test_obj._adjust_peak_bases()
    assert 'left_bases_adjusted' in test_obj.peaks['peak_props']
    assert 'right_bases_adjusted' in test_obj.peaks['peak_props']


def test_calculate_elution_volumes():
    test_obj.calculate_elution_volumes()
    assert 'elution_volume' in test_obj.peaks