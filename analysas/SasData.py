import numpy as np

class SasData:

    """

    Class for storing small-angle scattering data and performing simple reduction manipulations.

    Attributes:
    -----------

    q : array_like
        Units : Angstrom^-1
        One-dimensional scattering vector.
    Iq : array_like
        Units : cm^-1
        One-dimensional scattering intensity.
    dIq : array_like, optional
        Default value is 0.
        Units : cm^-1
        One-dimensional error associated with Iq.
    dq : array_like, optional
        Default value is 0.
        Units : Angstrom^-1
        One-dimensional error associated with q.
    label : string
        Sample name or label for the dataset.
    description: string, optional
        Extended description regarding the sample, dataset, data collection procedures, etc.

    """


    def __init__(self, label, description=''):

        """

        Constructor for SasData class.

        Parameters:
        -----------

        label : string
            Sample name or label for the dataset.

        Optional Parameters:
        --------------------

        description : string
            Extended description regarding the sample, dataset, data collection procedures, etc.

        """

        # confirm data inputs

        assert (type(label) is str), "Argument 'label' must be of type 'str', not " + type(label).__name__
        self.label = label

        assert (type(description) is str), "Argument 'description' must be of type 'str', not " + type(description).__name__
        self.description = description

        # initialize primary attributes

        self.q = np.array([])
        self.Iq = np.array([])
        self.dIq = np.array([])
        self.dq = np.array([])

        # initialize private variables

        self.__config_tags = np.array([])  # saved configurations as specified by the user
        self.__remove_tags = np.array([])  # all data points that have been removed by the user or various methods
        self.__trim_tags = np.array([])    # data pionts that have been trimmed by the user, subset of 'remove_tags'
        self.__smear_tags = np.array([])   # data points removed due to smearing, subset of 'remove_tags'
        self.__q_raw = np.array([])        # raw q data as added by the user
        self.__Iq_raw = np.array([])       # raw Iq data as added by the user
        self.__dIq_raw = np.array([])      # raw dIq data as added by the user
        self.__dq_raw = np.array([])       # raw dq data as added by the user
        self.__incoh_bkgd = 0              # currently subtracted incoherent background

    #########################
    # Private Class Methods #
    #########################

    def _sort_data(self):

        """ Sorts the data, ordered by the 'q' primary attribute."""

        sort_ind = np.argsort(self.q)
        self.q = self.q[sort_ind]
        self.Iq = self.Iq[sort_ind]
        self.dIq = self.dIq[sort_ind]
        self.dq = self.dq[sort_ind]

    def _get_filtered_data(self):

        """

        Returns a sorted copy of the raw data filtered by the remove_tags:
        (q_filtered, Iq_filtered, dIq_filtered, dq_filtered)

        """

        q_filtered = np.delete(self.__q_raw, self.__remove_tags)
        Iq_filtered = np.delete(self.__Iq_raw, self.__remove_tags)
        dIq_filtered = np.delete(self.__dIq_raw, self.__remove_tags)
        dq_filtered = np.delete(self.__dq_raw, self.__remove_tags)

        return (q_filtered, Iq_filtered, dIq_filtered, dq_filtered)

    def _update_remove_tags(self):

        """
        Updates the 'remove_tags' private variable by compiling unique values from:
        trim_tags
        smear_tags

        """

        self.__remove_tags = np.unique(np.concatenate((self.__trim_tags, self.__smear_tags)))

    ###################################################
    # Methods to Add Scattering Data in Various Forms #
    ###################################################

    def add_data(self, q, Iq, dIq=None, dq=None, config_tag=None):

        """

        Adds scattering data to the current instance of the SasData class.

        Parameters:
        -----------

        q : array_like
            Units : Angstrom^-1
            One-dimensional scattering vector.
        Iq : array_like
            Units : cm^-1
            One-dimensional scattering intensity.
            Array length should match that of q.

        Optional Parameters:
        --------------------

        dIq : array_like or float
            Default value is 0.
            Units : cm^-1
            One-dimensional array of error associated with Iq.
            Array length should match that of q (and Iq).
            A float can be passed to apply the same value for all Iq.
        dq : array_like, optional
            Default value is 0.
            Units : Angstrom^-1
            One-dimensional error associated with q.
            Array length should match that of q (and Iq).
            A float can be passed to apply the same value for all Iq.
        config_tag : string
            Label used to specify a configuration or specific subset of the data.
            This can be used to later call that specific subset of data.

        """

        # confirm data type, dimensions and length of argument 'q'
        try:
            q = np.array(q).astype(float)
        except:
            print("Argument 'q' must be array-like with float or int values.")
            raise
        assert (q.ndim == 1), "Argument 'q' must be one-dimensional array."
        q_length = len(q)

        # confirm data type, dimensions and length of argument 'Iq'
        try:
            Iq = np.array(Iq).astype(float)
        except:
            print("Argument 'Iq' must be array-like with float or int values.")
            raise
        assert (Iq.ndim == 1), "Argument 'Iq' must be one-dimensional array."
        assert (len(Iq) == q_length), "Argument 'Iq' must be of same length as 'q'."

        # confirm data type, dimensions and length of optional argument dIq
        if dIq is None:
            dIq = np.full((q_length), 0)
        try:
            dIq = np.array(dIq).astype(float)
        except:
            print("If provided, argument 'dIq' must be array-like with float or int values, or a single float value.")
            raise
        assert (dIq.ndim == 1), "If provided, argument 'dIq' must be one-dimensional array."
        if len(dIq) > 1:
            assert (len(dIq) == q_length), "If provided, argument 'dIq' must be of same length as 'q' (and 'Iq')."
        else:
            dIq = np.full(q_length, dIq[0])

        # confirm data type, dimensions and length of optional argument dq
        if dq is None:
            dq = np.full((q_length), 0)
        try:
            dq = np.array(dq).astype(float)
        except:
            print("If provided, argument 'dq' must be array-like with float or int values, or a single float value.")
            raise
        assert (dq.ndim == 1), "If provided, argument 'dq' must be one-dimensional array."
        if len(dq) > 1:
            assert (len(dq) == q_length), "If provided, argument 'dq' must be of same length as 'q' (and 'Iq')."
        else:
            dq = np.full(q_length, dq[0])

        # storing raw data in the respective private variables
        self.__q_raw = np.concatenate((self.__q_raw, q), axis = 0)
        self.__Iq_raw = np.concatenate((self.__Iq_raw, Iq), axis = 0)
        self.__dIq_raw = np.concatenate((self.__dIq_raw, dIq), axis = 0)
        self.__dq_raw = np.concatenate((self.__dq_raw, dq), axis = 0)

        # storing configuration information
        if config_tag is not None:
            assert (type(config_tag) is str), "If provided, argument 'config_tag' must be a string type."
        self.__config_tags = np.concatenate((self.__config_tags, np.full(q_length, config_tag)), axis=0)

        # adding the data to primary attributes
        # any current data treatments (e.g. remove incoherent background) will be applied
        self.q = np.concatenate((self.q, q), axis=0)
        self.Iq = np.concatenate((self.Iq, Iq - self.__incoh_bkgd), axis=0)
        self.dIq = np.concatenate((self.dIq, dIq), axis=0)
        self.dq = np.concatenate((self.dq, dq), axis=0)
        self._sort_data()

    def add_data_from_ABS(self, filepath, config_tag=None, skiprows=None):

        """

        Reads scattering data in ABS file format generated by the NIST Center for Neutron Research Igor reduction macros.

        This method looks for the following line within the file as an indication to the start of the data:
            'The 6 columns are | Q (1/A) | I(Q) (1/cm) | std. dev. I(Q) (1/cm) | sigmaQ | meanQ | ShadowFactor|'

        Parameters:
        -----------

        filepath : string
            Specifies the fielpath to the .ABS file.

        Optional Parameters:
        --------------------

        config_tag : string
            Label used to specify a configuration or specific subset of the data.
            This can be used to later call that specific subset of data.
        skiprows : int
            Default is None.
            Specifies how many rows to skip before the data begins.
            This is helpful with the file cannot be parsed easily.

        """

        if skiprows is None:
            num = 1
            file = open(filepath, 'r')
            for line in file:
                if line[:17] == 'The 6 columns are':
                    skiprows = num
                num += 1
            file.close()
            assert (type(skiprows) is int), "Method was unable to detect data in the ABS file at " + filepath

        assert (type(skiprows) is int), "An integer value must be provided for the 'skiprows' argument."
        try:
            data = np.loadtxt(filepath, skiprows = skiprows)
        except:
            print("Error extracting the data from the ABS file.")
            raise

        q = data[:,0]
        Iq = data[:,1]
        dIq = data[:,2]
        dq = data[:,3]

        self.add_data(q, Iq, dIq=dIq, dq=dq, config_tag=config_tag)

    def add_data_from_USANS(self, filepath, config_tag=None, skiprows=4):

        """

        Reads desmeared USANS data in txt file format generated by the NIST Center for Neutron Research Igor reduction macros.


        Parameters:
        -----------

        filepath : string
            Specifies the fielpath to the USANS file.

        Optional Parameters:
        --------------------

        config_tag : string
            Label used to specify a configuration or specific subset of the data.
            This can be used to later call that specific subset of data.
        skiprows : int
            Default value is 4.
            Specifies how many rows of header to skip at the beginning of the file.

        """
        assert (type(skiprows) is int), "An integer value must be provided for the 'skiprows' argument."

        try:
            data = np.loadtxt(filepath, skiprows=skiprows)
        except:
            print("Error extracting  data from the USANS file.")
            raise

        q = data[:,0]
        Iq = data[:,1]
        dIq = data[:,2]
        dq = None

        self.add_data(q, Iq, dIq=dIq, dq=dq, config_tag=config_tag)

    #################################
    # Methods for Miscellaneous Use #
    #################################

    def get_config_list(self):

        """

        Returns:
        --------

        config-list : list
            The list of configurations available to sort the data.
            Items in this list can be used as arguments for the get_config_data() method.

        """

        config_return = [x for x in self.__config_tags.tolist() if type(x) is not None]
        config_return = set(config_return) # get unique values
        config_return = list(config_return)
        return config_return


    def get_config_data(self, configs):

        """

        Enables access to subsets of the scattering data for a specific instrument and/or configuration.

        Parameters:
        -----------

        configs : str, iterable of str
            List of labels used to select a subset of the data based on instrument and/or configuration.

        Returns:
        --------

        (q, Iq, Iq_error, q_error) : tuple
            Data filtered by the provided configurations in 'configs'.

        """

        # confirming data inputs
        if type(configs) is str:
            configs = [configs]
        check_configs = self.get_config_list()
        for config in configs:
            assert (type(config) is str), "All configurations must be in string format."
            assert (config in check_configs), "Configuration '" + config + "' is not an available configuration."

        # filtering based on the user-specified configuration tags
        config_tags = np.delete(self.__config_tags, self.__remove_tags)
        filter_ind = np.isin(config_tags, configs)

        q, Iq, Iq_error, q_error = self._get_filtered_data()
        q_config = q[filter_ind]
        Iq_config = Iq[filter_ind] - self.__incoh_bkgd
        dIq_config = dIq_config[filter_ind]
        dq_config = dq_config[filter_ind]

        return (q_config, Iq_config, dIq_config, dq_config)

    def get_raw_data(self):

        """

        Returns a copy of the raw data:
        (q_raw, Iq_raw, dIq_raw, dq_raw)

        """

        q_raw = np.copy(self.__q_raw)
        sort_ind = np.argsort(q_raw)

        q_raw = q_raw[sort_ind]
        Iq_raw = np.copy(self.__Iq_raw)[sort_ind]
        dIq_raw = np.copy(self.__dIq_raw)[sort_ind]
        dq_raw = np.copy(self.__dq_raw)[sort_ind]

        return(q_raw, Iq_raw, dIq_raw, dq_raw)

    def reset_data(self):

        """
        Resets all data treatments (i.e. returns primary data to the raw data)

        """

        q, Iq, dIq, dq = self.get_raw_data()

        self.q = q
        self.Iq = Iq
        self.dIq = dIq
        self.dq = dq

        self.__incoh_bkgd = 0
        self.__smear_tags = np.array([])
        self.__trim_tags = np.array([])
        self._update_remove_tags()

    ###############################################
    # Methods to Remove the Incoherent Background #
    ###############################################

    def get_incoh_bkgd(self):

        """ Returns the currently subtracted incoherent background value. """

        return self.__incoh_bkgd

    def remove_inc_background(self, min_pts=20, set_val=None, threshold = 0.01):

        """

        Removes the flat, incoherent background at high-Q values from the scattering data.
        By calling this function, the incoherent background will be re-calculated and the primary attributes updated.

        Parameters:
        -----------

        min_pts : int, optional
            Default value is 20.
            This is the minimum number of points to average at the high-Q limit when determining the value of incoherent background.
        set_val : float, optional
            If specified, the user can overwrite the function's fitting procedure and specify the value of incoherent background to subtract.
            This can be useful if data is irregular at high-Q.
        threshold : float, optional
            Default value is 0.01.
            Increase this value to include more points at high-Q when calculating the incoherent background value.
            At high values, this can lead to an artificially high value of the incoherent background.

        """

        # assert that data input types are correct
        assert (type(min_pts) is int), "The argument 'min_pts' must be an integer value."
        try:
            threshold = float(threshold)
        except:
            print("The 'threshold' argument should be a numeric value only.")
            raise

        # refresh the primary attribute data by getting the filtered raw data
        self.q, self.Iq, self.dIq, self.dq = self._get_filtered_data()
        self._sort_data()

        # user can specify a background value if desired otherwise the data is fit at high-q
        if set_val is not None:
            try:
                set_val = float(set_val)
            except:
                print("The 'set_val' argument should be a numeric value only.")
                raise
            else:
                self.__incoh_bkgd = set_val

        else:
            avg, std = None, None
            for i in range(min_pts, len(self.Iq)):
                avg_new = np.average(self.Iq[-min_pts-i:])
                std_new = np.std(self.Iq[-min_pts-i:])
                if i == min_pts:
                    # the first average is always accepted
                    avg = avg_new
                    std = std_new
                else:
                    # only add a point and accept the new average if it is within the threshold
                    if np.absolute(avg_new - avg) > threshold*std:
                        break
                    else:
                        avg = avg_new
                        std = std_new
            self.__incoh_bkgd = avg

        self.Iq = self.Iq - self.__incoh_bkgd

    def reset_incoh_bkgd(self):

        """ Resets the incoherent background removal (i.e. sets incoherent background to 0). """

        self.__incoh_bkgd = 0
        self.q, self.Iq, self.dIq, self.dq = self._get_filtered_data()
        self._sort_data()

    ##################################
    # Methods to Remove Smeared Data #
    ##################################

    def remove_smeared_data(self, configs=[], threshold = 0.75):

        """

        Removes data points that have a q_error/q value greater than the threshold (0.75 by default).
        Note that this will not reset any previous removals of smeared data.

        To reset (clear) all smeared data removals, utilize the 'reset_smeared_data' method.


        Optional Parameters:
        -----------

        configs : iterable
            Specify which configurations smeared data should be removed from.
            Default behavior will be to remove smeared data from all configurations.

        threshold : float
            Default is 0.75.
            Any data points with dq/q greater than the threshold will be removed.

        """

        if type(configs) is str:
            configs = [configs]

        if len(configs) == 0:
            configs = self.get_config_list()
        else:
            check_configs = self.get_config_list()
            if len(check_configs) == 0:
                print ("There are no configurations associated with the data for: " + self.label)
                raise
            else:
                for config in configs:
                    assert (type(config) is str), "All configurations must be in string format."
                    assert (config in check_configs), "Configuration '" + config + "' is not an available configuration."

        try:
            threshold = float(threshold)
        except:
            print("The 'threshold' argument must be a numeric value.")
            raise

        dq_q = np.absolute(self.__dq_raw/self.__q_raw)

        smear_filter = np.where(dq_q > threshold, True, False)
        config_filter = np.isin(self.__config_tags, configs)
        all_filters = np.logical_and(smear_filter, config_filter)
        new_smear_tags = np.arange(0,len(dq_q))[all_filters]

        self.__smear_tags = np.unique(np.concatenate((new_smear_tags, self.__smear_tags), axis=0))
        self._update_remove_tags()
        self.q, self.Iq, self.dIq, self.dq = self._get_filtered_data()
        self.Iq = self.Iq - self.__incoh_bkgd
        self._sort_data()

    def reset_smeared_data(self):

        """ Resets all smeared data removals. """

        self.__smear_tags = np.array([])
        self._update_remove_tags()
        self.q, self.Iq, self.dIq, self.dq = self._get_filtered_data()
        self.Iq = self.Iq - self.__incoh_bkgd
        self._sort_data()

    ########################
    # Methods to Trim Data #
    ########################

    def trim_q_range(self, q_min, q_max):

        """

        Removes data points whose scattering vector is outside the range [q_min, q_max]


        Parameters:
        -----------

        q_min : float
            Lower edge of acceptable scattering vector range, inclusive.
        q_max : float
            Upper edge of acceptable scattering vector range, inclusive.
            q_max must be greater than q_min.

        """

        try:
            q_min = float(q_min)
            q_max = float(q_max)
        except:
            print("Scattering vector range [q_min, q_max] must be provided as numeric values only")

        assert (q_min < q_max), "Ensure q_min is less than q_max."

        low_trim = np.where(self.__q_raw < q_min)[0]
        high_trim = np.where(self.__q_raw > q_max)[0]
        full_trim = np.concatenate((low_trim, high_trim), axis=0)

        self.__trim_tags = np.concatenate((self.__trim_tags, full_trim), axis=0)

        self._update_remove_tags()
        self.q, self.Iq, self.dIq, self.dq = self._get_filtered_data()
        self.Iq = self.Iq - self.__incoh_bkgd
        self._sort_data()

    def trim_points(self, trim_points):

        """

        Removes data points specified by index.
        Indexing starts at 0 and increases with scattering vector q.

        Indices must be provided in a list, but an item in the list can
        include a tuple specifying the range (inclusive, exclusive).

        Example:
            trim_points = [0, 5, (10, 14), 20]
            The points removed would have the indices:
            [0, 5, 10, 11, 12, 13, 20]

        Parameters:
        -----------

        trim_points : list
            Data indices to remove.

        """

        format_trim = []
        for ind in trim_points:
            if type(ind) is int:
                format_trim.append(ind)
            elif type(ind) is tuple:
                format_trim += list(range(ind[0],ind[1]))
            else:
                pass

        update_trim = []
        for ind in format_trim:
            q_loc = np.where(self.__q_raw == self.q[ind])[0]
            Iq_loc = np.where(self.__Iq_raw == self.Iq[ind] + self.get_incoh_bkgd())[0]
            loc = np.intersect1d(q_loc, Iq_loc)
            for val in loc:
                update_trim.append(val)

        full_trim = np.unique(np.array(update_trim))
        self.__trim_tags = np.concatenate((self.__trim_tags, full_trim), axis=0)

        self._update_remove_tags()
        self.q, self.Iq, self.dIq, self.dq = self._get_filtered_data()
        self.Iq = self.Iq - self.__incoh_bkgd
        self._sort_data()

    def trim_config(self, configs):

        """

        Removes data points specified by configuration.

        Parameters:
        -----------

        configs : str, iterable of str
            List of labels used to remove a subset of the data based on instrument and/or configuration.

        """

        if type(configs) is str:
            configs = [configs]

        check_configs = self.get_config_list()
        for config in configs:
            assert (type(config) is str), "All configurations must be in string format."
            assert (config in check_configs), "Configuration '" + config + "' is not an available configuration."

        filter_ind = np.isin(self.__config_tags, configs)
        full_trim = np.arange(0, len(self.__q_raw))[filter_ind]

        self.__trim_tags = np.concatenate((self.__trim_tags, full_trim), axis=0)

        self._update_remove_tags()
        self.q, self.Iq, self.dIq, self.dq = self._get_filtered_data()
        self.Iq = self.Iq - self.__incoh_bkgd
        self._sort_data()

    def trim_points(self, trim_points):

        """

        Removes data points specified by index.
        Indexing starts at 0 and increases with scattering vector q.

        Indices must be provided in a list, but an item in the list can
        include a tuple specifying the range (inclusive, exclusive).

        Example:
            trim_points = [0, 5, (10, 14), 20]
            The points removed would have the indices:
            [0, 5, 10, 11, 12, 13, 20]

        Parameters:
        -----------

        trim_points : list
            Data indices to remove.

        """

        format_trim = []
        for ind in trim_points:
            if type(ind) is int:
                format_trim.append(ind)
            elif type(ind) is tuple:
                format_trim += list(range(ind[0],ind[1]))
            else:
                pass

        update_trim = []
        for ind in format_trim:
            q_loc = np.where(self.__q_raw == self.q[ind])[0]
            Iq_loc = np.where(self.__Iq_raw == self.Iq[ind] + self.get_incoh_bkgd())[0]
            loc = np.intersect1d(q_loc, Iq_loc)
            for val in loc:
                update_trim.append(val)

        full_trim = np.unique(np.array(update_trim))
        self.__trim_tags = np.concatenate((self.__trim_tags, full_trim), axis=0)

        self._update_remove_tags()
        self.q, self.Iq, self.dIq, self.dq = self._get_filtered_data()
        self.Iq = self.Iq - self.__incoh_bkgd
        self._sort_data()

    def reset_trim(self):

        """
        Removes any data trimming that has been performed.
        """

        self.__trim_tags = np.array([])

        self._update_remove_tags()
        self.q, self.Iq, self.dIq, self.dq = self._get_filtered_data()
        self.Iq = self.Iq - self.__incoh_bkgd
        self._sort_data()
