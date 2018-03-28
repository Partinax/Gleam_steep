from astropy.io import fits
from optparse import OptionParser
from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u
from astroquery.skyview import SkyView
from astropy.table import Table, Column
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
import os
import matplotlib.pyplot as plt

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

def format_filter_table(results_table): # formats the query return table to be ready for writing to fis table, also removs unwanted object types such as star
    temp_otype = results_table.field('OTYPE')
    results_table['MAIN_ID'] = results_table['MAIN_ID'].astype(str) # convert these table columns dtypes from object to str so that they can be saved to a .fits
    results_table['OTYPE'] = results_table['OTYPE'].astype(str)
    table_row_num = 0
    columns_removed = 0
    for row in temp_otype:
        if row == "Star":
            results_table.remove_row(table_row_num - columns_removed)
            columns_removed += 1
        table_row_num += 1
    return results_table

def write_to_fits(table, name): # writes each individual cross match instance to an individual .fits table
    table.write(name + ".fits", overwrite=True)
    table.pprint()

def plot_data_list(fits, matches, cat, object_id, name): # want to also plot the cross match table
    f1 = cat.field('int_flux_088')
    f2 = cat.field('int_flux_118')
    f3 = cat.field('int_flux_154')
    f4 = cat.field('int_flux_200')
    f5 = cat.field('alpha')

    y = [f1[object_id], f2[object_id], f3[object_id], f4[object_id]]
    x = [88,118,154,200]
    fig = plt.figure()
    sed = fig.add_subplot(231)
    sed.plot(x, y)
    sed.text(118, f2[object_id], r'$\alpha=$' + str(f5[object_id]))

    matches.write('table', format='latex') # purpose of all this mess is to convert the Table object to a latex format, then read that back in as a single string
    with open('table', 'r') as myfile:
        text = myfile.read().replace("\\begin{table}", '') # various conversions to get matplot lib to read as a table
        text = text.replace( "\end{table}", '')
        text = "r'''"+text+"'''"
    table = fig.add_subplot(232)
    table.text(0, 0, text, size=12)
    index = 3
    for hdul in fits: # need a good way to determine if the image is NVSS TGSS or DSS
        wcs = WCS(hdul[0].header)
        image = fig.add_subplot(2, 3, index)
        fig.add_subplot(projection=wcs)
        image.imshow(hdul[0].data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
        image.grid(color='white', ls='solid')
        image.set_xlabel('Right Ascension J2000') # can't add labels directly to the subplot object?
        image.set_ylabel('Declination J2000')
        index += 1
    fig.savefig(name+".png")


def get_images_list(RA, DEC):
    pos = RA + " " + DEC
    image_list=[]
    if DEC > -53.0:
        tgss = SkyView.get_images(position=pos, survey='TGSS ADR1')
        image_list.append(tgss[0])
    if DEC > -40.0:
        nvss = SkyView.get_images(position=pos, survey='NVSS')
        image_list.append(nvss[0])

    dss = SkyView.get_images(position=pos, survey='DSS')
    image_list.append(dss[0])
    return image_list

def get_images_save(RA, DEC):
    pos = RA + " " + DEC
    if DEC > -53.0:
        tgss = SkyView.get_images(positio=pos, survey='TGSS ADR1')
        tgss.writeto(pos+"_TGSS.fits")
    if DEC > -40.0:
        nvss = SkyView.get_images(position=pos, survey='NVSS')
        nvss.writeto(pos + "_NVSS.fits")
    dss = SkyView.get_images(position=pos, survey='DSS')
    dss.writeto(pos + "_DSS.fits")



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
catalogue='filtered_red_selection_fin.fits'
t = Table.read(catalogue, format = 'fits')


ob_row_num = 0
for i in cat_data:
    tempra = str(i[0])
    tempdec = str(i[2])
    string2 = tempra + " " + tempdec
    string = tempra + "d " + tempdec + "d"
    c = coordinates.SkyCoord(string, frame='icrs')
    r = 2 * u.arcminute
    results_table = customSimbad.query_region(c, radius=r)
    fits = get_images_list(tempra, tempdec)
    #image_list = SkyView.get_images(position=string2, survey=['NVSS', 'TGSS ADR1', 'DSS'])
    if type(results_table) is Table:
        results_table = format_filter_table(results_table)
        temp_z = results_table.field('Z_VALUE')
        temp_asec = results_table.field('ANG_DIST')
        seperation_dist = calc_sep_distance(temp_z, temp_asec)
        sep_dist=Column(name="SEP_DIST(Mpc)", data=seperation_dist)
        results_table.add_column(sep_dist)
        write_to_fits(results_table, string2)
        #print image_list
        write_to_txt(results_table, string2, options.output)
        plot_data_list(fits, results_table, t, ob_row_num, string2)
    ob_row_num += 1



