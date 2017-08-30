# data_cleaning_for_disp

The code was used to clean up active ingredient names in the Drug_Product_Ingredients table from [Canada Vigilance Adverse Reaction Online Database](https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database/canada-vigilance-adverse-reaction-online-database-data-structures.html)

The goal was to recognize active ingredients in salt forms, recognize synonyms for an active ingredient and correct typos,in order to map
active ingredient names to ATC code, NHP names or vaccine products.

The ultimate goal was to improve accuracy in the calculation of disproportionality analysis (the new_cv_drug_rxn.R file is the primer for
later calculation ie.master_table_pt,master_table_hlt creation for the Disproportionality App)
