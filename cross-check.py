from astropy.io import fits
from optparse import OptionParser
from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u
from astroquery.skyview import SkyView
from astropy.table import Table, Column

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


def write_to_txt(results_table, string2, output):
    f = open(output + ".txt", 'a')
    f.write('####################################### ' + string2 + ' #######################################' + '\n')
    f.write(str(results_table))
    f.write('\n')
    f.close()

def format_filter_table(results_table):
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

def write_to_fits(table, name):
    table.write(name + ".fits", overwrite=True)
    table.pprint()


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

ob_row_num = 0
for i in cat_data:
    tempra = str(i[0])
    tempdec = str(i[2])
    string2 = tempra + " " + tempdec
    string = tempra + "d " + tempdec + "d"
    c = coordinates.SkyCoord(string, frame='icrs')
    r = 2 * u.arcminute
    results_table = customSimbad.query_region(c, radius=r)
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
    ob_row_num += 1



