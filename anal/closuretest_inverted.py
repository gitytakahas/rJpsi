from ROOT import TFile, Double
import math

procs = ['bc_others', 'bc_jpsi_dst', 'bc_jpsi_tau', 'bg_ul']


for year in ['2016', '2017', '2018']:
    
    mc = {}
    data = {}

    for reg in ['sr', 'sb']:

        file = TFile(year + '_' + reg + '_None/datacard/tau_rhomass_unrolled_var.root')
        
#        print file 

        # data - Bc MC 

        tot = 0.
        err = 0.
        check = 0.

        for proc in procs:
            hist = file.Get(reg + '/' + proc)

            err_ = Double(-1.)
            val = hist.IntegralAndError(1, hist.GetXaxis().GetNbins(), err_)

            factor = 1.
            if proc =='bc_jpsi_tau':
                factor = 0.25

    
            if proc=='bg_ul':
                tot += val*factor
            else:
                tot -= val*factor
                check += val*factor


            err += math.pow(err_*factor, 2)


        hist_data = file.Get(reg + '/data_obs')

        err_data = Double(-1.)
        val_data = hist_data.IntegralAndError(1, hist_data.GetXaxis().GetNbins(), err_data)

        mc[reg] = {'val':tot, 'err':math.sqrt(err)}
        data[reg] = {'val':val_data, 'err':err_data}

        
        print reg, year, check/val_data #reg, 'pred:', tot, '+/-', math.sqrt(err), 'obs:', val_data, '+/-', err_data, 'diff=', float(tot - val_data)/val_data

    ratio_data = data['sr']['val']/data['sb']['val']
    ratio_mc = mc['sr']['val']/mc['sb']['val']
    
    err_data = math.sqrt(math.pow(data['sr']['err']/data['sr']['val'],2) + math.pow(data['sb']['err']/data['sb']['val'],2))*ratio_data
    err_mc = math.sqrt(math.pow(mc['sr']['err']/mc['sr']['val'],2) + math.pow(mc['sb']['err']/mc['sb']['val'],2))*ratio_mc

    print year, 'data = {0:.4f}'.format(ratio_data), '+/- {0:.4f}'.format(err_data)
    print year, 'mc   = {0:.4f}'.format(ratio_mc), '+/- {0:.4f}'.format(err_mc)
    print year, 'diff   = {0:.4f}'.format((ratio_data - ratio_mc)/ratio_mc)
    
