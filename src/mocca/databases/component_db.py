#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:32:36 2021

@author: haascp
"""
from mocca.databases.base import BaseDatabase

class ComponentDatabase(BaseDatabase):
    def __init__(self):
        self.components = []

    def _average_peak_spectrum(self, peak):
        """
        Calculates mean spectrum over peak from left to right border.
        """
        return np.average(peak.dataset.data[:, peak.left:peak.right+1],
                          axis=0).tolist()

    def _average_spectra_over_peaks(self, peaks):
        """
        Calculates mean spectrum of a list of peaks with averaged spectrum.
        """
        spectra_list = []
        for peak in peaks:
            peak_spectrum = self._average_peak_spectrum(peak)
            spectra_list.append(peak_spectrum)
        return np.average(np.array(spectra_list), axis=0)

    def _average_peak_retention_times(self, peaks):
        """
        Calculates mean retention indices of a list of peaks.
        """
        left = int(round(sum([peak.left for peak in peaks]) / len(peaks))),
        right = int(round(sum([peak.right for peak in peaks]) / len(peaks))),
        maximum = int(round(sum([peak.maximum for peak in peaks]) / len(peaks)))
        return left, right, maximum

    def _get_valid_peaks(self, peaks):
        """
        Filters list of peaks for pure and unsaturated peaks with a compound_id.
        """
        return [peak for peak in peaks if (peak.pure is True and
                                           peak.saturation is False and
                                           peak.compound_id is not None)]

    # TODO: implement different filter options for relevant_peaks (e.g. by date)
    def _filter_peaks(self, peaks, condition):
        if condition is None:
            return peaks
        else:
            return peaks

    def _sort_peaks_by_compound(self, peaks):
        """
        Creates dict with unique compound_id as keys and a list of corresponding
        peaks as values.
        """
        compound_dict = {}
        for peak in peaks:
            if peak.compound_id not in compound_dict:
                compound_dict[peak.compound_id] = []
            compound_dict[peak.compound_id].append(peak)
        return compound_dict

    def _create_component_from_peaks(self, compound_id, peaks):
        """
        Creates a component object based on the given peaks
        """
        mean_left, mean_right, mean_maximum = self._average_peak_retention_times(peaks)
        mean_spectrum = self._average_spectra_over_peaks(peaks)
        return Component(compoound_id=compound_id,
                         left=mean_left,
                         right=mean_right,
                         maximum=mean_maximum,
                         spectrum=mean_spectrum,
                         created_from=peaks
                         )

    def update(self, peak_database, condition=None):
        """
        Creates components from the given peak database. Optionally, a condition
        can be given to filter peaks.
        """
        # clear database to fill it with components
        self._delete_all_components()

        # only peaks with assigned compound_id which are not saturated and pure
        valid_peaks = self._get_valid_peaks(peak_database.peaks)

        # only peaks satisfying filter condition
        filtered_peaks = self._filter_peaks(valid_peaks, condition)

        # collect peaks by compound
        compound_dict = self._sort_peaks_by_compound(filtered_peaks)

        # create components out of compound dict
        for compound_id, peaks in compound_dict.items():
            self.components.append(self._create_component_from_peaks(compound_id,
                                                                     peaks))

    # for addition of unknown compounds in reaction runs if initialization runs
    # are not available anymore
    def insert_by_compound_id(self, peak_database, compound_id, condition=None):
        """
        Inserts component in existing component list. If component with given
        compound_id already exists, it will be overwritten.
        """
        if compound_id in self.components:
            self.components = [comp for comp in self.components if
                               comp.compound_id is not compound_id]
            logging.debug("Component with compound_id {} already existed in "
                          "components. New component owerwrites old component."
                          "".format(compound_id))

        compound_peaks = [peak for peak in peak_database.peaks if
                          peak.compound_id == compound_id]

        # only peaks with assigned compound_id which are not saturated and pure
        valid_peaks = self._get_valid_peaks(compound_peaks)

        # only peaks satisfying filter condition
        filtered_peaks = self._filter_peaks(valid_peaks, condition)

        self.components.append(self._create_component_from_peaks(compound_id,
                                                                 filtered_peaks))

    def __contains__(self, compound_id):
        """Allows statements such as "'Component 1' in component_database"
        to see if that Component is inside the database"""
        return compound_id in [comp.compound_id for comp in self.components]

    def __getitem__(self, compound_id):
        """Allows statements such as "component_database['Component 1']" to access
           that Component"""
        if compound_id not in self:
            return next(comp for comp in self.components if
                        comp.compound_id == compound_id)
        else:
            raise AttributeError("{} not found in Component Database!"
                                 "".format(compound_id))

    def __iter__(self):
        """Yields all components inside the database.
           Allows statements such as "for component in component_database:" """
        for component in self.components:
            yield component