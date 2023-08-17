import pandas as pd

"""
Testing filter 1 function

bc_str_df: dataframe with barcode, STR, count
"""


def BarcodeToCGenotypeProcessor(bc_str_df):
    """
    One barcode might map to multiple similar
	(or different) STR sequences. This function
	performs genotypeing similar to HipSTR to obtain
	a best guess STR sequence for each barcode

    Arguments
    --------
    bc_str_df: pandas Dataframe
    Returns
    --------

    """
    unique_bc = set(bc_str_df["BC"])
    for bc in unique_bc:
        BarcodeToGenotye(bc)
    return None


bc_str_df = pd.DataFrame({
    "BC": ["GGTTGGTGTTTGTCTTGTTT"]* 6 + ['CTCGGTTGTATATTCTATAG'] * 2,
    "STR": ['chr12:114684488-114684503_Human-STR-339270_AGCG_12_AT_0',
            'chr12:114684488-114684503_Human-STR-339270_AGCG_13_AT_0',
             'chr12:106957456-106957468_Human-STR-333288_A_9_A_0',
              'chr12:114684488-114684503_Human-STR-339270_AGCG_11_AT_0',
               'chr12:114684488-114684503_Human-STR-339270_AGCG_14_AT_0',
                'chr12:114684488-114684503_Human-STR-339270_AGCG_10_AT_0'] +
                ['Human_STR_228160_p5', 'Human_STR_228160_ref'],
    "Count": [231, 20, 2, 17, 4, 1, 13, 1]

})

print (bc_str_df)

bc_str_cleaned = BarcodeToCGenotype(bc_str_df)
print (bc_str_cleaned)
