"""
Additional Features:
    - Peak Deconvolution (for better peak integration, etc.)

Misc. Thoughts:
    - Try / Except is nice, but leads to black box smoothing. Probably better to throw error?
    - Would want to think longer about data structuring in the obj, esp processed_data.
    - How does data smoothing impact accuracy of peak integration w.r.t. to peak height?
    - How to handle 'step' values, with the one n.a.?
    - How to handle data consistency, since we only support one set of smoothed data / peaks at a time?
      Could break that into another object, but need to consider the workflow.
"""

import matplotlib.pyplot as plt
from mpl_interactions import panhandler, zoom_factory
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import numpy as np
from BaselineRemoval import BaselineRemoval


class LCData:
    def __init__(self, metadata, data):
        """Initializes an LCData object.

        Args:
            metadata (dict): Dictionary of the metadata from the LC run.
            data (list[float, str, float]): Array containing the raw data as [time, step, value].
        """

        self.metadata = metadata
        self.data = data

    # --- Convenient getters
    def get_time(self):
        """Getter for time values.

        Returns:
            list[float]: 1-D array of the time values.
        """

        return [row[0] for row in self.data]

    def get_step(self):
        """Getter for step values.

        Returns:
            list[str]: 1-D array of the step values as strings.
        """

        return [row[1] for row in self.data]

    def get_value(self):
        """Getter for signal values.

        Returns:
            list[float]: 1-D array of the signal values.
        """

        return [row[2] for row in self.data]

    # --- Data processing
    def process_raw_chromatogram(self, sigma=1.0):
        """Perfroms data smoothing and baseline correction using gaussian 1D smoothing algorithm and Zhang fit algorithm,
        respectively. Assigns results to object property, ydata_processed.

        Scipy Docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html
        BaselineRemoval Docs: https://pypi.org/project/BaselineRemoval/

        Args:
            sigma (float, optional): Standard deviation for Gaussian kernel. Defaults to 1.
        """

        ydata = self.get_value()
        ydata_smoothed = gaussian_filter1d(ydata, sigma)
        base_obj = BaselineRemoval(ydata_smoothed)
        ydata_baseline_corrected = base_obj.ZhangFit()
        self.ydata_processed = ydata_baseline_corrected

    # --- Peak methods
    def detect_peaks(self, width=0.0, height=0.0, distance=3):
        """Detects peaks using the processed ydata for the LC run. If no processed data is present, it will be generated
        using default values.

        Args:
            width (float, optional): Required width of peaks in samples. Defaults to 0.
            height (float, optional): Required height of peaks. Defaults to 0.
            distance (int, optional): Required minimal horizontal distance (>= 1) in samples between neighbouring peaks.
              Smaller peaks are removed first until the condition is fulfilled for all remaining peaks. Defaults to 3.

        In Scipy, width, height, and distance can all be None, but the corresponding metadata will not be present, 
        hence the numeric defaults here.

        Scipy Docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
        Threshold/Prominence: https://stackoverflow.com/questions/76786670/what-is-the-difference-between-threshold-and-prominence-in-scipy-find-peaks
        left/right_bases often double counts sections of the chromatogram, as noted here: https://github.com/scipy/scipy/issues/19232.
        Hence, we generated adjusted peak_bases which ensures no overlap between peaks.
        """

        try:
            ydata = self.ydata_processed
        except AttributeError:
            self.process_raw_chromatogram(sigma=1)

        peaks, peak_props = find_peaks(
            ydata, width=width, height=height, distance=distance
        )
        self.peaks = {"peak_index": peaks, "peak_props": peak_props}
        self._adjust_peak_bases()

    def integrate_peaks(self):
        """Integrates the detected peaks as an approximation of elution volume.

        This will use left/right_bases_adjusted data to avoid double counting peak area.
        Peak integration can be improved by adding peak deconvolution to object.
        """

        xdata = self.get_time()
        try:
            ydata = self.ydata_processed
            peak_props = self.peaks["peak_props"]
        except AttributeError:
            self.detect_peaks()
            ydata = self.ydata_processed
            peak_props = self.peaks["peak_props"]

        lower = peak_props["left_bases_adjusted"]
        upper = peak_props["right_bases_adjusted"]

        peak_areas = [
            np.trapz(ydata[lb:ub], x=xdata[lb:ub]) for lb, ub in zip(lower, upper)
        ]
        self.peaks["peak_areas"] = peak_areas

    def _adjust_peak_bases(self):
        """Function for adjusting peak boundaries to avoid overlapping peak sections.

        Crunchy method that does not yield perfect results. For non-overlapping sections of peaks, a 'peak' will be 
           added for the baseline inbetween sections, but since this is used primarily for peak integration and these
           peaks will have ~0 area, it is okay for a first pass.
        Should be fixed / updated for production work.
        """

        lb = self.peaks.get("peak_props").get("left_bases")
        rb = self.peaks.get("peak_props").get("right_bases")
        peak_bases_sorted = sorted(list(set(np.concatenate((lb, rb)))))
        lb_new = peak_bases_sorted[0:-1]
        rb_new = peak_bases_sorted[1::]

        self.peaks["peak_props"]["left_bases_adjusted"] = lb_new
        self.peaks["peak_props"]["right_bases_adjusted"] = rb_new

    def calculate_elution_volumes(self):
        """Calculates the elution volumes corresponding to each peak.

        Assuming all intejected fluid exits within the timeframe of the experiment, than the area of a peak is 
        proportional to the elution volume (also assuming linear signal intensity, independent of species).
        """

        injection_volume = float(
            self.metadata["Injection Information"]["Injection Volume (µL)"]
        )
        peak_areas = self.peaks["peak_areas"]
        total_area = np.sum(peak_areas)
        elution_volumes = [a / total_area * injection_volume for a in peak_areas]
        self.peaks["elution_volume"] = elution_volumes


def parse_lc_textfile(file):
    """Parses a Liquid Chromatography txt data file (matching the example data format) into an LCData object.

    Args:
        file (str): Path to text file.

    Returns:
        LCData: Object containing chromatogram data and metadata.

    Note for file format encoding: https://stackoverflow.com/questions/17912307/u-ufeff-in-python-string
    """

    f = open(file, "r", encoding="utf-8-sig")
    lines = f.readlines()

    section = None
    metadata = {}
    data = []

    for i, line in enumerate(lines):
        line_split = line.replace("\n", "").split("\t")

        if line == "\n":
            continue

        # Determines sections based on line length, and adds metadata to appropriate section.
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

        # Cromatogram Data Section is used to generate a data list, not a metadata dictionary.
        if section == "Chromatogram Data":
            data.append(line_split)

    # Removes units from chromatogram data section and adds to metadata
    data_units = data[0]
    metadata["Chromatogram Data Information"]["Units"] = data_units

    # Fixes datatypes for chromatogram data
    data_new = []
    for row in data[1::]:
        row_new = [float(row[0]), row[1], float(row[2])]
        data_new.append(row_new)

    return LCData(metadata, data_new)


def plot_raw_vs_processed_data(obj, xlim=None, ylim=None):
    """Generates a plot of the raw data vs the processed data (smoothing and baseline correction).

    Args:
        obj (LCData): LCData object containing data and metadata.
        xlim (None|list[float, float], optional): Sets the xlim for the plot. Defaults to None.
        ylim (None|list[float, float], optional): Sets the ylim for the plot. Defaults to None.
    """

    xdata = obj.get_time()
    ydata_raw = obj.get_value()
    try:
        ydata_processed = obj.ydata_processed
    except AttributeError:
        obj.process_raw_chromatogram()
        ydata_processed = obj.ydata_processed

    # Enable scroll to zoom with the help of MPL
    # Interactions library function like ioff and zoom_factory.
    plt.clf()
    with plt.ioff():
        figure, axis = plt.subplots()

    plt.xlabel(obj.metadata["Chromatogram Data Information"]["Units"][0])
    plt.ylabel(obj.metadata["Chromatogram Data Information"]["Units"][2])
    plt.plot(xdata, ydata_raw, label="Raw Data")
    plt.plot(xdata, ydata_processed, label="Processed Data")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.legend()

    # Enable scrolling and panning with the help of MPL
    # Interactions library function like panhandler.
    disconnect_zoom = zoom_factory(axis)
    pan_handler = panhandler(figure)
    display(figure.canvas)


def plot_peak_data(obj, xlim=None, ylim=None, peak_lines=True, peak_areas=False):
    """Generates a plot of the processed data and denotes the peak positions.

    Args:
        obj (LCData): LCData object containing data and metadata.
        xlim (None|list[float, float], optional): Sets the xlim for the plot. Defaults to None.
        ylim (None|list[float, float], optional): Sets the ylim for the plot. Defaults to None.
        peak_lines (bool, optional): Indicates if lines showing peaks will be drawn. Defaults to True.
        peak_areas (bool, optional): Indicates if inegrated peak areas will be drawn. Defaults to False.
    """

    try:
        peaks = obj.peaks["peak_index"]
        peak_props = obj.peaks["peak_props"]
    except AttributeError:
        obj.detect_peaks()
        peaks = obj.peaks["peak_index"]
        peak_props = obj.peaks["peak_props"]

    xdata = obj.get_time()
    ydata = obj.ydata_processed

    # Enable scroll to zoom with the help of MPL
    # Interactions library function like ioff and zoom_factory.
    plt.clf()
    with plt.ioff():
        figure, axis = plt.subplots()

    plt.xlabel(obj.metadata["Chromatogram Data Information"]["Units"][0])
    plt.ylabel(obj.metadata["Chromatogram Data Information"]["Units"][2])
    plt.plot(xdata, ydata, "-", label="Processed Data")
    if peak_lines:
        plt.vlines(
            [xdata[i] for i in peaks],
            [0 for _ in peaks],
            [h for h in peak_props["peak_heights"]],
            colors="red",
            linestyle="-.",
        )
    if peak_areas:
        lower = obj.peaks["peak_props"]["left_bases_adjusted"]
        upper = obj.peaks["peak_props"]["right_bases_adjusted"]
        for i, (lb, ub) in enumerate(zip(lower, upper)):
            if i % 2 == 0:
                color = "lightgreen"
            else:
                color = "darkgreen"
            plt.fill_between(xdata[lb:ub], ydata[lb:ub], color=color)

    plt.xlim(xlim)
    plt.ylim(ylim)

    # Enable scrolling and panning with the help of MPL
    # Interactions library function like panhandler.
    disconnect_zoom = zoom_factory(axis)
    pan_handler = panhandler(figure)
    display(figure.canvas)
