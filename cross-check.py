from astropy.io import fits
from optparse import OptionParser
from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u
from astroquery.skyview import SkyView
from astropy.table import Table, Column
from astropy.wcs import WCS
import astropy.wcs as wcs
from astropy.utils.data import get_pkg_data_filename
import os
import matplotlib as mpl
mpl.use('Agg')  # So does not use display
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import numpy as np
import gleam_client as GC
import matplotlib.patches as patches
import montage_wrapper as montage
import math
from astroquery.vizier import Vizier


def generate_powerlaw(start, stop, C1, C2):
    Y = []
    X = []
    x_range = np.arange(start, stop)
    print(stop)
    print(x_range[len(x_range)-1])
    for i in range(0, len(x_range)):
        Y.append(C1 + C2*x_range[i])
        X.append(x_range[i])

    return X,Y



def calc_sep_distance(z, asec):
    sd = []
    H = 73.45  # uncertainty of 1.66 Riess18
    for i in range(0,len(z)):
        if z[i] > 0.0:
            v = z[i] * 3.*10.**5.
            d = v / H
            rad = asec[i] * 4.848*10**(-6.) # double check this conversion
            sd.append(rad*d)
        else:
            sd.append("--")
    return sd


def write_to_txt(results_table, string2, output): #writes all the collected cross match data to a single text file
    f = open(output + ".txt", 'a')
    f.write('####################################### ' + string2 + ' #######################################' + '\n')
    f.write(str(results_table))
    f.write('\n')
    f.close()
    #results_table.pprint()

def format_filter_table(results_table): # formats the query return table to be ready for writing to fis table, also removs unwanted object types such as star
    temp_otype = results_table.field('OTYPE')
    results_table['MAIN_ID'] = results_table['MAIN_ID'].astype(str) # convert these table columns dtypes from object to str so that they can be saved to a .fits
    results_table['OTYPE'] = results_table['OTYPE'].astype(str)
    #hresults_table['COORDINATES'] = results_table['COORDINATES'].astype(str)
    table_row_num = 0
    columns_removed = 0
    for row in temp_otype:
        if row == "Star":
            results_table.remove_row(table_row_num - columns_removed)
            columns_removed += 1
        table_row_num += 1
    temp_otype = results_table.field('OTYPE')
    if len(temp_otype)>10:
        results_table.remove_rows(slice(10, len(temp_otype)-1,1))
    return results_table

def truncate_table(table, length):
    if len(table)>length:
        table.remove_rows

def check_matches(results_table):
    temp_otype = results_table.field('OTYPE')
    #hresults_table['COORDINATES'] = results_table['COORDINATES'].astype(str)
    table_row_num = 0
    columns_removed = 0
    match = False
    i = 0
    while match == False and i<len(temp_otype):

        if temp_otype[i] == "Galaxy":
            match = True
        if temp_otype[i] == "GinCl":
            match == True
        if temp_otype[i] == "BClG":
            match == True

        i += 1
    return match

def write_to_fits(table, name): # writes each individual cross match instance to an individual .fits table
    table.write(name + ".fits", overwrite=True)


def plot_data_list(fits_list, matches, cat, object_id, object_name, flux_results, galactic, alpha): # Handles all the final fomratting and data presentation for the table output
    #### preping metadata and table format ####

    gleam_RA = cat.field('ra_088')[object_id]
    gleam_DEC = cat.field('dec_088')[object_id]
    f1 = cat.field('int_flux_088')
    f2 = cat.field('int_flux_118')
    f3 = cat.field('int_flux_154')
    f4 = cat.field('int_flux_200')
    f5 = cat.field('alpha')
    err_f1 = cat.field('err_int_flux_088')
    err_f2 = cat.field('err_int_flux_118')
    err_f3 = cat.field('err_int_flux_154')
    err_f4 = cat.field('err_int_flux_200')
    peak = cat.field('peak_flux_154')[object_id]

    first = {"marker": "+", "linestyle": "None", "color": "red"}
    second = {"marker": "o", "linestyle": "None", "color": "red"}
    third = {"marker": "x", "linestyle": "None", "color": "red"}
    fourth = {"marker": "*", "linestyle": "None", "color": "red"}
    fifth = {"marker": "v", "linestyle": "None", "color": "red"}
    gleam_marker = {"marker": ".", "linestyle": "None", "color": "blue"}

    if type(matches) is Table:
        ra = matches.field('RA')
        dec = matches.field('DEC')
        matches.remove_column('MAIN_ID')
        matches.remove_column('ANG_DIST')
        matches.remove_column('Z_VALUE')
        matches.remove_column('RA')
        matches.remove_column('DEC')
        marker = ['']*(len(ra))
        marker_list = ['+', 'o', 'x', '*', 'v']
        for j in range(len(ra)):
            if type(marker[j]) is not None and j < 5:
                marker[j] = marker_list[j]
        marker = Column(name="MARKER", data=marker)
        results_table.add_column(marker)
        matches.pprint()
        matches.write('table', format='latex') # purpose of all this mess is to convert the Table object to a latex format, then read that back in as a single string
        with open('table', 'r') as myfile:

            matches_list = myfile.readlines()#.replace("\\begin{table}","")  # various convmatches.remove_column('RA')ersions to get matplot lib to read as a table

            del matches_list[0]
            del matches_list[len(matches_list) - 1]
            matches_string = ''
            for line in matches_list:
                matches_string += line


        matches_string = '%r' %matches_string
        matches_string = matches_string.replace('\\n', '')
        matches_string = matches_string.replace('\\\\', '\\')
        print matches_string

    else:
        matches_string = "none"
        ra = []

    icons = [first, second, third, fourth, fifth]



    #### plotting data ####
    fig = plt.figure()

    ## plot SED ##

    y = np.array([f1[object_id], f2[object_id], f3[object_id], f4[object_id]])
    x = np.array([88, 118, 154, 200])
    x_clean, y_clean = clean_list(x, y)

    log_y = np.log(y)
    log_x = np.log(x)


    if gleam_DEC > -40.0:
        freq = np.log(1400)
        key = "S1.4"
    else:
        freq = np.log(847)
        key = "St"
    fit = np.polyfit(np.log(x_clean), np.log(y_clean), 1)
    fit_x,fit_y = generate_powerlaw((np.log(70)),(8), fit[1], fit[0])
    print(fit)
    sed = fig.add_subplot(231)
    err = np.array([err_f1[object_id], err_f2[object_id], err_f3[object_id], err_f4[object_id]])
    np.log(err)
    sed.errorbar(log_x, log_y, yerr= err,fmt="b.")
    sed.plot(fit_x, fit_y)
    try:
        flux = np.log(flux_results[0][key]/1000)
        print(object_id)
        print(flux_results[0][key])
        sed.plot(freq, flux, marker="X")
    except:
        flux = None
    sed.text(np.log(88), np.log(f1[object_id]), r'$\alpha=$' + str(f5[object_id])+r'$\alpha=$' + str(alpha))
    #sed.set_yscale('log')
    #sed.set_xscale('log')
    sed.set_xlabel('LOG(Frequency (MHz))')
    sed.set_ylabel('LOG(Flux (Jy))')


    ## plot matches table ##
    mpl.rc('text', usetex=True)
    table = fig.add_subplot(232)
    table.text(0, 0, matches_string, size=10)

    ## plot TGSS, NVSS/SUMSS and DSS images ##
    index = 3 # will always start at the third subplot position
    for i in range(0, len(fits_list), 2):
        name = fits_list[i]
        hdul = fits_list[i+1]
        wcs = WCS(hdul[0].header)
        image = fig.add_subplot(2, 3, index, projection=wcs)
        c = SkyCoord(str(gleam_RA) + " " + str(gleam_DEC), unit=(u.deg, u.deg))
        x, y = wcs.wcs_world2pix(c.ra.degree, c.dec.degree, 0)
        image.plot(x, y, **gleam_marker)
        for j in range(len(ra)): # adds co-ordinates markers to images for 5 closest simbad matches
            if type(ra[j]) is not None and j < 5: # checks to see if markers haven't gone over and ra is long enough for each marker
                c = SkyCoord(str(ra[j])+" "+str(dec[j]), unit=(u.hourangle, u.deg))

                x, y = wcs.wcs_world2pix(c.ra.degree, c.dec.degree, 0) # how does this even work? It doesn't know the FOV of the images
                image.plot(x, y, **icons[j])

        image.grid(color='white', ls='solid')
        image.set_xlabel('Right Ascension J2000')
        image.set_ylabel('Declination J2000')
        if name == 'TGSS':
            image_temp = hdul[0].data
            vmin = 1 * np.std(image_temp)
            vmax = 4 * vmin
            image.imshow(hdul[0].data, vmin=-vmax, vmax=vmax, origin='lower') #dynamically applies image scaling
            image.set_title("TGSS")

        if name == 'NVSS':
            image_temp = hdul[0].data
            vmin = 1 * np.std(image_temp)
            vmax = 4 * vmin
            image.imshow(hdul[0].data, vmin=-vmax, vmax=vmax, origin='lower')
            image.set_title("NVSS")

        if name == 'DSS':
            image_temp = hdul[0].data
            vmin = 1 * np.std(image_temp)
            vmax = 10 * vmin
            image.imshow(image_temp, vmin=vmin, vmax=vmax, origin='lower')
            image.set_title("DSS")

        if name == 'SUMSS':
            image_temp = hdul[0].data
            vmin = 1 * np.std(image_temp)
            vmax = 4 * vmin
            image.imshow(image_temp, vmin=-vmin, vmax=vmax, origin='lower')
            image.set_title("SUMSS")

        index += 1

    ## plot gleam image ##
    ra = gleam_RA
    dec = gleam_DEC
    ang_size = 0.50
    freq_low = ['072-103', '103-134', '139-170']
    projection = 'SIN'
    dl_dir = '/home/matt/project/Data/temp/'

    GC.vo_get(ra, dec, ang_size, proj_opt=projection, freq=freq_low, download_dir=dl_dir)
    gleam_name_r = dl_dir + GC.create_filename(ra, dec, ang_size, freq_low[0])
    gleam_name_g = dl_dir + GC.create_filename(ra, dec, ang_size, freq_low[1])
    gleam_name_b = dl_dir + GC.create_filename(ra, dec, ang_size, freq_low[2])


    #fits.writeto("red_test.fits", red_data)
    #fits.writeto("green_test.fits", green_data)
    try:
        montage.mGetHdr(gleam_name_b, "temp.txt")
        blue_hdul = fits.open(gleam_name_b)
        blue_data = blue_hdul[0].data

        new_gleam_name_r = dl_dir + "_repr_" + GC.create_filename(ra, dec, ang_size, freq_low[0])
        montage.reproject(gleam_name_r, new_gleam_name_r, header="temp.txt", exact_size=True)
        red_hdul = fits.open(new_gleam_name_r)
        red_data = red_hdul[0].data

        new_gleam_name_g = dl_dir + "_repr_" + GC.create_filename(ra, dec, ang_size, freq_low[1])
        montage.reproject(gleam_name_g, new_gleam_name_g, header="temp.txt", exact_size=True)
        green_hdul = fits.open(new_gleam_name_g)
        green_data = green_hdul[0].data
        RGB = np.dstack([red_data, green_data, blue_data])
        RGB -= np.nanmin(RGB)# 0 to 1 scale
        RGB /= (np.nanmax(RGB) + 3*peak)/2.
        wcs = WCS(blue_hdul[0].header)

    except: # if there is some data lacking to create the RGB image, then the program just defaults to a red image (80 MHz).
        red_hdul = fits.open(gleam_name_r)
        red_data = red_hdul[0].data
        RGB = red_data
        wcs = WCS(red_hdul[0].header)


    image = fig.add_subplot(2, 3, index, projection=wcs)
    image.imshow(RGB, origin = 'lower')
    image.grid(color='white', ls='solid')
    image.set_xlabel('Right Ascension J2000')
    image.set_ylabel('Declination J2000')
    image.set_title("GLEAM")
    c = SkyCoord(str(gleam_RA) + " " + str(gleam_DEC), unit=(u.deg, u.deg))
    x, y = wcs.wcs_world2pix(c.ra.degree, c.dec.degree, 0)
    image.add_patch(patches.Rectangle((12.5, 12.5), 25, 25, fill=False, edgecolor="white"))
    image.plot(x, y, **gleam_marker)



    ## save image to png ##
    fig.set_size_inches((13, 8.5), forward=False)
    if galactic == False:
        plt.savefig("/home/matt/project/Data/aggregates/No_matches/"+str(object_id+1) +', '+ object_name+".png", dpi=500)
    else:
        plt.savefig("/home/matt/project/Data/aggregates/Matches/" + str(object_id + 1) + ', ' + object_name + ".png", dpi=500)
    plt.close()


def get_images_list(RA, DEC, radius):
    pos = RA + " " + DEC
    image_list=[]
    DEC = float(DEC)
    if DEC > -53.0:
        tgss = SkyView.get_images(position=pos, survey='TGSS ADR1', radius=radius)
        image_list.append('TGSS')
        image_list.append(tgss[0])
    if DEC > -40.0:
        nvss = SkyView.get_images(position=pos, survey='NVSS', radius=radius)
        image_list.append('NVSS')
        image_list.append(nvss[0])
    else:
        sumss = SkyView.get_images(position=pos, survey='SUMSS 843 MHz', radius=radius)
        image_list.append('SUMSS')
        image_list.append(sumss[0])

    dss = SkyView.get_images(position=pos, survey='DSS', radius=radius)
    image_list.append('DSS')
    image_list.append(dss[0])
    return image_list

def get_images_save(RA, DEC, id): # if exists use this method instead of downloading new one.
    pos = RA + " " + DEC
    DEC = float(DEC)
    if DEC > -53.0:
        tgss = SkyView.get_images(position=pos, survey='TGSS ADR1')
        tgss[0].writeto('/home/matt/project/Data/temp/'+ str(id) +"_" +pos+"_TGSS.fits", overwrite = True)
    if DEC > -40.0:
        nvss = SkyView.get_images(position=pos, survey='NVSS')
        nvss[0].writeto('/home/matt/project/Data/temp/'+ str(id) +"_" + pos + "_NVSS.fits", overwrite = True)
    else:
        sumss = SkyView.get_images(position=pos, survey='SUMSS 843 MHz')
        sumss[0].writeto('/home/matt/project/Data/temp/'+ str(id) +"_" + pos + "_SUMSS.fits", overwrite = True)
    dss = SkyView.get_images(position=pos, survey='DSS')
    dss[0].writeto('/home/matt/project/Data/temp/'+ str(id) +"_" + pos + "_DSS.fits", overwrite = True)

def clean_list(x, y):
    removed = 0
    for i in range(0, len(y)):
        if (math.isnan(y[i-removed])):
            y = np.delete(y, i - removed)
            x = np.delete(x, i - removed)
            removed += 1
    return x, y

def calc_alpha(x, y):
    sums = 0
    alpha = 0
    y =y/1000.
    x =x*1000000.
    if len(x) > 1:
        for i in range(1, len(x)):
            alpha_i = (np.log(y[i])-np.log(y[i-1]))/(np.log(x[i])-np.log(x[i-1]))
            sums += alpha_i
            alpha = sums/(len(x)-1)
            #print alpha_i

    return alpha



################# __MAIN__ ######################
usage = "Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('--catalogue', type="string", dest="catalogue",
                  help="The filename of the catalogue you want to read in.", default=None)
parser.add_option('--output', type="string", dest="output",
                  help="The name of the output catalogue.", default=None)
(options, args) = parser.parse_args()



cat = fits.open(options.catalogue)
cat_data = cat[1].data
customSimbad = Simbad()
customSimbad.add_votable_fields('dist(asec)')
customSimbad.add_votable_fields('otype')
customSimbad.add_votable_fields('z_value')
customSimbad.remove_votable_fields('coordinates')
customSimbad.add_votable_fields('ra')
customSimbad.add_votable_fields('dec')
catalogue='filtered_red_selection_fin.fits'
t = Table.read(catalogue, format = 'fits')
end = len(cat_data)

ob_row_num = 31 # starting position in data table
for i in range(ob_row_num, end): # this counter is a bit redundant currently, should update
    try:

        alpha_list = []
        f1 = t.field('int_flux_088')
        f2 = t.field('int_flux_118')
        f3 = t.field('int_flux_154')
        f4 = t.field('int_flux_200')
        y = np.array([f1[i], f2[i], f3[i], f4[i]])
        x = np.array([88, 118, 154, 200])
        x, y = clean_list(x, y)
        alpha = calc_alpha(x, y)
        alpha_list = alpha_list.append(alpha)
        tempra = str(t.field('ra_088')[i])
        tempdec = str(t.field('dec_088')[i])
        string2 = tempra + " " + tempdec
        string = tempra + "d " + tempdec + "d"
        c = coordinates.SkyCoord(string, frame='icrs')
        r = 4 * u.arcminute # simbad query radius.
        image_radius = 15 * u.arcminute
        #v = Vizier(columns=['phot.flux*', 'em_radio*'])
        flux_results = Vizier.query_region(c,radius=1.0 * u.arcminute, catalog=['NVSS','SUMSS'])
        results_table = customSimbad.query_region(c, radius=r)
        fits_list = get_images_list(tempra, tempdec, image_radius)
        get_images_save(tempra, tempdec, i + 1)
        #image_list = SkyView.get_images(position=string2, survey=['NVSS', 'TGSS ADR1', 'DSS'])
        galactic = False
        if type(results_table) is Table:
            results_table = format_filter_table(results_table)
            galactic = check_matches(results_table)
            temp_z = results_table.field('Z_VALUE')
            temp_asec = results_table.field('ANG_DIST')
            seperation_dist = calc_sep_distance(temp_z, temp_asec)
            sep_dist=Column(name="SEP DIST(Mpc)", data=seperation_dist)
            results_table.add_column(sep_dist)
            #write_to_fits(results_table, string2)
            #print image_list
            write_to_txt(results_table, string2, options.output)
        plot_data_list(fits_list, results_table, t, i, string2, flux_results, galactic, alpha)
    except:
        print("Failed Object "+str(i+1))
        write_to_txt(str(i+1),"failed object","failed.txt")

#np.savetxt('alpha_list.txt', np.c_[alpha_list])
histogram = np.histogram(alpha_list)
plt.plot(histogram[1], histogram[0])


