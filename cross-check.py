from astropy.io import fits
from optparse import OptionParser
from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u
from astroquery.skyview import SkyView
from astropy.table import Table, Column

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
    image_list = SkyView.get_image_list(position=string2, survey=['NVSS', 'TGSS ADR1', 'DSS'])
    if type(results_table) is Table:
        temp_otype = results_table.field('OTYPE')
        results_table['MAIN_ID']=results_table['MAIN_ID'].astype(str) # convert these table columns dtypes from object to str so that they can be saved to a .fits
        results_table['OTYPE'] = results_table['OTYPE'].astype(str)
        table_row_num = 0
        columns_removed = 0
        for x in temp_otype:
            if (x == "Star"):
                results_table.remove_row(table_row_num - columns_removed)
                columns_removed += 1
            table_row_num += 1
        name = str(i[0]) + '_' + str(i[2])
        results_table.pprint()
        results_table.write('new.fits')
    print image_list
    ob_row_num += 1
