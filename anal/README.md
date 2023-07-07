# NOTE
# Step 1. to 4 is NOT necessary unless you change the dataset
# Step 5. is to produce plots locally 
# Step 6. is to produce plots in batch (note that --min option is specified in job_template.sh. Remove this option if you want to produce all plots ... but it takes long time ...)



###########################################################
# 1. produce inclusive plots with minimum option and no weighting applied for inclusive J/psi bkg.
# after enabling 'inclusive' category in draw.py, run,

# normal selections
> python draw.py --year 2018 --min && python draw.py --year 2017 --min && python draw.py --year 2016 --min

# inverted tau selections
> python draw.py --year 2018 --min --inv && python draw.py --year 2017 --min --inv && python draw.py --year 2016 --min --inv



###########################################################
# 2. calculate normalization factor for the inclusive J/psi Backgrounds

> python calculateNorm4inclusiveJpsiBkg.py 

# multiply printed-out additional factors to the variable ``bkg_data_sf'' in draw.py 


###########################################################
# 3. run again draw.py with updated scale factor for inslucive J/psi Bkg. You might want to remove --min option to produce all plots in inclusive category. 

# normal selections
> python draw.py --year 2018 --min && python draw.py --year 2017 --min && python draw.py --year 2016 --min

# inverted tau selections
> python draw.py --year 2018 --min --inv && python draw.py --year 2017 --min --inv && python draw.py --year 2016 --min --inv

# To check, you can do again,

> python calculateNorm4inclusiveJpsiBkg.py 

# to see all the SFs to the inclusive J/psi BG is now close to 1.


###########################################################
# 4. derive correction factors for the BDT by fitting [0, 3.5] range for 2017/18 and [0, 3.] for 2016.

> python norm_BDT.py


###########################################################
# 5. ploduce all plots locally. Disable 'inclusive' category and enable other categories (typically, sr, sb, lp and gap). Then, run draw.py with bkg corrections applied and all plots options. 

# normal selections
> python draw.py --year 2018 -w && python draw.py --year 2017 -w && python draw.py --year 2016 -w

# inverted tau selections
> python draw.py --year 2018 --inv -w && python draw.py --year 2017 --inv -w && python draw.py --year 2016 --inv -w


###########################################################
# 6. produce plots in batch. 
# Disable 'inclusive' category and enable other categories (typically, sr, sb, lp and gap). Then, run draw.py with bkg corrections applied and ``--min'' plots options. 

> python getDataset_simultaneous.py


###########################################################
# 7. produce datacards

> python createFinalDatacard_simultaneousfit.py
> python createFinalDatacard_simultaneousfit.py --scale (needed for the combine studies)


###########################################################
# 8. produce validation plots

> python shape_gap_sb.py 
> python shape_comparison_display.py 
> python createFinalDatacard.py
> python createFinalDatacard.py --inv


###########################################################

# 9. update plots to AN (change target directory inside the script!).
Don't forget to update your repository before doing so!

> sh copy2ANrepo.sh
