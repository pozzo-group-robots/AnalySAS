from .SasData import *

def subtract_sas(A, B, label = None, description = '', scale_B = 1, threshold = 0.05):

    """
    Subtracts the intensity of two scattering profiles:
    I(Q)_A - scale_B * I(Q)_B = I(Q)_C

    It is recommended that data to be subtracted follows identical reduction protocols.
    This function does not do a fit of the 'B' dataset, and so it is only able to subtract
    the scattering intensities if the scattering vector 'q' is within a threshold value.

    The threshold value is calculated as:
    | q_B - q_A | / q_A

    Parameters:
    -----------

    A : SasData object
    B : SasData object

    Optional Parameters:
    --------------------

    label : string
        The sample label that will get applied to the returned SasData object C.
        Default behavior will utilize the label of SasData object A.

    description : string
        The extended description that will get applied to the returned SasData object C.
        Regardless of whether this is specified, default behavior will be to append:
        "A subtraction of [A.label] - [scale_B]*[B.label] was performed."

    scale_B : float
        Default value is 1.
        The scale factor can be used to adjust the contribution of scattering from B.
        This is especially important when subtracting solvent contribution from a mixture,
        and the solute concentration is high enough where volume displacement applies.

    threshold : float
        Data is subtracted only for q values from two datasets are within the threshold.

    Returns:
    --------

    C : SasData object
        Subtracted scattering data.

    """

    try:
        scale_B = float(scale_B)
    except:
        print ("The 'scale_B' argument must be a numeric value.")
        raise

    try:
        threshold = float(threshold)
    except:
        print ("The 'threshold' argument must be a numeric value.")


    sub_ind = []
    for q in A.q:
        diffs = np.absolute(B.q - q)
        sub_ind.append(np.argmin(diffs))

    Iq_C = A.Iq - scale_B * B.Iq[sub_ind]
    dIq_C = np.sqrt(A.dIq**2 + (scale_B * B.dIq[sub_ind])**2)

    q_diffs = np.absolute(A.q - B.q[sub_ind])/A.q
    filter_diffs = np.where(q_diffs <= threshold)[0]

    q_C = A.q[filter_diffs]
    Iq_C = Iq_C[filter_diffs]
    dIq_C = dIq_C[filter_diffs]
    dq_C = A.dq[filter_diffs]

    filter_negs = np.where(Iq_C >= 0)[0]

    q_C = q_C[filter_negs]
    Iq_C = Iq_C[filter_negs]
    dIq_C = dIq_C[filter_negs]
    dq_C = dq_C[filter_negs]

    if label is None:
        label = A.label
    else:
        try:
            label = str(label)
        except:
            print("WARNING: Your provided label could not be converted to a string, instead using: " + A.label)
            label = A.label

    assert (type(description) is str), "The provided description must be a string."
    if (len(description) > 0 and description[-1] != ' '):
        description += ' '
    description += ("NOTE: A subtraction of '" + A.label +"' - " + str(scale_B) + " * '" + B.label + "' was performed.")

    C = SasData(label, description = description)
    C.add_data(q_C, Iq_C, dIq = dIq_C, dq = dq_C)

    return C
