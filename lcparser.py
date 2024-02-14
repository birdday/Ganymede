"""
TODO:
    - Check naming conventions for classes, methods, and functions.
    - Doctstrings
    - Tests
    - Setup.py
"""

import matplotlib.pyplot as plt
from mpl_interactions import ioff, panhandler, zoom_factory
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import numpy as np

def parse_lc_textfile(file):
    f = open(file, "r",  encoding="utf-8-sig")
    # https://stackoverflow.com/questions/17912307/u-ufeff-in-python-string
    lines = f.readlines()

    section = None
    metadata = {}
    data = []

    for i, line in enumerate(lines):
        line_split = line.replace("\n", "").split("\t")

        if line == "\n":
            continue

        if len(line_split) == 1:
            section = line_split[0].replace(":", "")
            if section != "Chromatogram Data":
                metadata[section] = {}
            continue
        
        if not section:
            key, val = line_split[0], line_split[1]
            metadata[key] = val

        if section and section != "Chromatogram Data":
            key, val = line_split[0], line_split[1]
            metadata[section][key] = val

        if section == "Chromatogram Data":
            data.append(line_split)

    data_units = data[0]
    metadata["Chromatogram Data Information"]["Units"] = data_units

    data_new = []
    for row in data[1::]:
        row_new = [float(row[0]), row[1], float(row[2])]
        data_new.append(row_new)

    return LCData(metadata, data_new)


class LCData:
    def __init__(self, metadata, data):
        self.metadata = metadata
        self.data = data

    def get_time(self):
        return [row[0] for row in self.data]
    
    def get_step(self):
        return [row[1] for row in self.data]
    
    def get_value(self):
        return [row[2] for row in self.data]
    
    def detect_peaks(self, width=0, height=0):
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
        # https://stackoverflow.com/questions/76786670/what-is-the-difference-between-threshold-and-prominence-in-scipy-find-peaks
        ydata = self.get_value()
        # width and height default to None in scipy, but need a value for necessary metadata to be returned.
        peaks, peak_props = find_peaks(ydata, width=width, height=height)
        self.peaks = peaks
        self.peak_props = peak_props
    
    def integrate_peaks(self):
        xdata = self.get_time()
        ydata = self.get_value()
        peak_indicies = self.peaks
        peak_props = self.peak_props

        lower = peak_props["left_bases"]
        upper = peak_props["right_bases"]

        peak_areas = [np.trapz(ydata[lb:ub], x=xdata[lb:ub]) for lb, ub in zip(lower, upper)]
        self.peak_areas = peak_areas

    def smooth_data(self, sigma=1):
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html
        ydata = self.get_value()
        self.values_smoothed = gaussian_filter1d(ydata, sigma)


def plot_raw_data(obj):
    xdata = obj.get_time()
    ydata = obj.get_value()

    # Enable scroll to zoom with the help of MPL
    # Interactions library function like ioff and zoom_factory.
    plt.clf()
    with plt.ioff():
        figure, axis = plt.subplots()

    plt.xlabel(obj.metadata["Chromatogram Data Information"]["Units"][0])
    plt.ylabel(obj.metadata["Chromatogram Data Information"]["Units"][2])
    plt.plot(xdata, ydata)

    # Enable scrolling and panning with the help of MPL
    # Interactions library function like panhandler.
    disconnect_zoom = zoom_factory(axis)
    pan_handler = panhandler(figure)
    display(figure.canvas)


def plot_peak_data(obj):
    xdata = obj.get_time()
    ydata = obj.get_value()
    try:
        peaks = obj.peaks
        peak_props = obj.peak_props
    except AttributeError:
        obj.detect_peaks()
        peaks = obj.peaks
        peak_props = obj.peak_props

    # Enable scroll to zoom with the help of MPL
    # Interactions library function like ioff and zoom_factory.
    plt.clf()
    with plt.ioff():
        figure, axis = plt.subplots()

    plt.xlabel(obj.metadata["Chromatogram Data Information"]["Units"][0])
    plt.ylabel(obj.metadata["Chromatogram Data Information"]["Units"][2])
    plt.plot(xdata, ydata)
    plt.vlines(
        [xdata[i] for i in peaks], [0 for _ in peaks],
        [h for h in peak_props["peak_heights"]],
        colors='red',
        linestyle='-.'
    )

    # Enable scrolling and panning with the help of MPL
    # Interactions library function like panhandler.
    disconnect_zoom = zoom_factory(axis)
    pan_handler = panhandler(figure)
    display(figure.canvas)
