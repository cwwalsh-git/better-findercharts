#!/usr/bin/env python
# coding: utf-8

# In[1]:


import logging
import pkg_resources
import astropy
from astropy import coordinates
from astropy.coordinates import SkyCoord, FK5
from astropy.wcs import wcs
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from astropy.table import Table
from astropy.visualization import PercentileInterval, AsinhStretch
import numpy as np
import numpy.ma as ma
import astroquery
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astroquery.skyview import SkyView
import matplotlib.pyplot as plt
import os
import csv
import urllib.request
from urllib.error import HTTPError
import time
import argparse
import requests
from io import StringIO
import itertools
import pandas as pd
from sklearn.neighbors import KDTree, NearestNeighbors
from sklearn.preprocessing import normalize
import textwrap
import warnings

def log_setup():
    """sets up logging file and std logging"""
    logging.getLogger("astroquery").setLevel(logging.WARNING)
    logger = logging.getLogger("MAIN")
    fname = os.getcwd() + "/finderchart.log"
    LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
    logging.basicConfig(filename=fname, level=logging.INFO, format=LOG_FORMAT, filemode='w+')
    stdlog = logging.StreamHandler()
    stdlog.setLevel(logging.INFO)
    stdformat = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    stdlog.setFormatter(stdformat)
    logging.getLogger('').addHandler(stdlog)
    return logger, stdlog
    
def parser():
    """reads arguments from command line"""
    parser = argparse.ArgumentParser(prog='testProgram', description='Returns finderchart image(s) from specificed locations, with overlaid data and markers.', epilog="Main surveys: DSS, SDSSg, Pan-STARRS1/PS1. \n")
    parser.add_argument('-p', '--pos', nargs="*", metavar='', help='SIMBAD-resolvable object name/coordinate position, or comma-separated list of positions.')
    parser.add_argument('-s', '--survey', nargs="*", default='DSS', help='Survey to retrieve images from, or list of surveys. -h to view all possible surveys.')
    parser.add_argument('-w', '--width', default=5, type=float, metavar='', help='Width and height of image in arcminutes')
    parser.add_argument('-m', '--no_markers', action='store_false', help='Boolean to show or hide overlaid markers on image.')
    parser.add_argument('-f', '--file_entry', default=0, type=int, metavar='', help='Default pos argument instead takes .csv file path as string. Each row in .csv should contain a SIMBAD-resolvable name or position. Set=2 if .csv has separate ra and dec columns, otherwise set=1.')
    parser.add_argument('-l', '--list_surveys', action='store_true', help='List all surveys instead of running program.')
    parser.add_argument('-q', '--quiet', help='Returns quiet output', action='store_true')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-a', '--use_all_surveys', action='store_true', help='Returns images from ALL surveys in survey list.')
    group.add_argument('-r', '--survey_redundant', action='store_true', help='If survey in list fails, will search again with next survey in list.')
    args = parser.parse_args()
    if args.list_surveys:
        print({'Optical: Pan-STARSS1 class=': ['ps1', "Pan-STARRS1"]})
        print(SkyView.survey_dict)
        raise SystemExit
    return args

def fileEntryParse(fname, fileEntry):
    """Allows for file entry of positions from .csv file"""
    logger.debug(f"{fname}")
    try:
        pos_df = pd.read_csv(fname, dtype=str)
    except FileNotFoundError as e:
        if ".csv" not in fname:
            logger.error(f".csv was not found in file name '{fname}'. Retrying by appending .csv to filename...")
            fname+=".csv"
        else:
            logger.error(f"No .csv file found at {fname}. Please move .csv into same directory as .py file.")
            raise e
        pos_df = pd.read_csv(fname, dtype=str)
    pos_df.columns = [chr(i) for i in range(ord('a'), ord('a')+pos_df.shape[1])]
    if fileEntry==1:
        positions = pos_df['a']
    elif fileEntry==2:
        positions = [*zip(pos_df['a'], pos_df['b'])]
        positions = [', '.join(i) for i in positions]
    else:
        raise ValueError(f"Value of fileEntry argument '{fileEntry}' is invalid. Please use False (list of position(s)), 1 (1-column csv entry), or 2 (2-column csv entry).")
        logger.error(f"Value of fileEntry argument '{fileEntry}' is invalid. Please use False (list of position(s)), 1 (1-column csv entry), or 2 (2-column csv entry).")
    return positions

def makeGaiaRegion(pos, img, wcsf, side):
    """get image file and retrieve rectangular coordinate region for Gaia query"""
    pos = pos.replace('(', '').replace(')', '')
    if ',' not in pos:
        pos = pos.replace(' ', ', ')
    shape = img.shape
    p1 = wcsf.pixel_to_world(0, 0)
    p2 = wcsf.pixel_to_world(shape[0], shape[0])
    p1, p2 = p1.to_string().split(' '), p2.to_string().split(' ')
    d1 = float(p1[0])-float(p2[0])
    d2 = float(p1[1])-float(p2[1])
    d = np.sqrt(d1**2 + d2**2)
    r = np.mean([float(p1[0]), float(p2[0])])
    dec = np.mean([float(p1[1]), float(p2[1])])
    center = f"{str(r)}, {str(dec)}"
    region = f"{center}, {d}, {d}"
    logger.info(f"Created search region for Gaia at ({region}).")
    return center, region

def regionGaiaQuery(pos, region):
    """find all objects in specified region within Gaia eDR3 catalog"""
    gaia_columns = "source_id, ra, dec, pmra, pmdec, ipd_frac_multi_peak"
    gaia_query_region = f"""SELECT
    {gaia_columns}, DISTANCE(POINT('ICRS', {str(pos)}), POINT('ICRS', ra, dec)) AS dist FROM gaiaedr3.gaia_source
    WHERE 1=CONTAINS(
        POINT('ICRS', ra, dec),
        BOX('ICRS', {region}))
    ORDER BY dist ASC"""
    try:
        gaia_job = Gaia.launch_job_async(gaia_query_region)
        results = gaia_job.get_results()
    except Exception as e:
        er = str(e).replace("\n", " ")
        logger.error(f"Gaia search failed with error '{er}'.")
        raise ConnectionError(f"No Gaia data found in search region. Gaia returned: '{er}'")
    gaiafname = 'gaia_results.csv'
    results.write(gaiafname, overwrite=True)  
    results_df = pd.read_csv(gaiafname)
    logger.info(f"Found {len(results)} source(s) from Gaia eDR3 in region.")
    return results, results_df

def epochMath(array, epoch):
    """applies epoch transform to coordinates (see epochAdjust)"""
    array = array.filled(fill_value=0)
    return np.array([*map(lambda x, epoch=epoch:(epoch-2016)*x/3_600_000, array)])

def epochAdjust(gaia_results, header, survey):
    """shift coordinates of each source by its proper motion to account for epoch difference between img and Gaia data"""
    gRes = gaia_results
    epochstr=''
    for value in header.values(): 
        try:
            if "epoch" in value.casefold():
                epochstr = ''.join(i for i in value if i.isnumeric())
        except AttributeError:
            continue
    if not epochstr:
        logger.info('Unable to get epoch data from header of image. Markers may be slightly misaligned with image.')
        return gRes
    epochs = textwrap.wrap(epochstr, 4)
    epochs = [*map(float, epochs)]
    epoch = np.mean(np.array(epochs))
    pmra_t = epochMath(gRes['pmra'], epoch)
    pmdec_t = epochMath(gRes['pmdec'], epoch)
    gRes['ra']+=pmra_t
    gRes['dec']+=pmdec_t
    return gRes
    
def plotPointandPM(gRes):
    """plot markers and proper motion arrows on all positions containing sources in Gaia eDR3 catalogue"""
    icrs = ax.get_transform('icrs')
    plt.plot(gRes['ra'], gRes['dec'], transform=icrs, linestyle='none', marker='o', markersize=15, markeredgecolor='white', markerfacecolor='none')
    for r, d, pmr, pmd in zip(gRes['ra'], gRes['dec'], gRes['pmra'], gRes['pmdec']):
        if pmr:
            ax.arrow(r, d, pmr/1440, pmd/1440, transform=icrs, head_width=0.0002, head_length=0.0002, color='w', width=0.00002)

def plotBinaries_fast(gRes):
    """connect likely binary pairs based on PM using dataframe/nearest neighbor optimization"""
    gResc = gRes[np.abs(gRes['pmra']) > 0]
    gResc = gResc.reset_index()
    gResPM = gResc.loc[:, ["pmra", "pmdec"]]
    try:
        gResPM_n = normalize(gResPM, norm='l2')
        nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree').fit(gResPM_n)
        distances, indices = nbrs.kneighbors(gResPM_n)
    except ValueError:
        return 1
    join = pd.DataFrame(np.concatenate((indices, distances), axis=1))
    join.columns = ['in1', 'in2', 'ign', 'dist']
    join = join[join.dist < 0.001]
    for in1, in2 in zip(join['in1'], join['in2']):
        a = np.array((gResc['ra'][in1], gResc['dec'][in1]))
        b = np.array((gResc['ra'][in2], gResc['dec'][in2]))
        if 0 < np.linalg.norm(a-b) < 0.05:
            plt.plot([gResc['ra'][in1], gResc['ra'][in2]], [gResc['dec'][in1], gResc['dec'][in2]], 
                     color='red', linestyle='dashed', linewidth=0.8, transform=ax.get_transform('icrs'))

def getImage(pos, survey, side, **kwargs):
    """use Skyview function to retrieve finderchart images"""
    if 'help' in survey.casefold():
        SkyView.list_surveys()
        raise("Hope that helped.")
    h = side*u.arcminute
    try:
        path = SkyView.get_images(position=pos, survey=survey, coordinates='ICRS', radius=h)
    except HTTPError as e:
        logger.error(e)
        logger.error(f"{survey} does not contain image at ({pos}).")
        raise ConnectionError(f"{survey} does not contain image at ({pos}).")
    fit = path[0]
    img = fit[0].data
    header = fit[0].header
    wcsf = wcs.WCS(fit[0])
    data = dict(header)
    pos = data["CRVAL1"], data["CRVAL2"]
    pos = str(pos)
    transform = AsinhStretch() + PercentileInterval(99.75)
    img = transform(img)
    logger.info(f"Retrieved image from {survey}.")
    return pos, img, header, wcsf

def PS1_geturl(ra, dec, arcmin_size=240, filters="i", output_size=None):
    """use PAN-STARRS1 image cutout service to get .fits url for finderchart"""
    size = int(arcmin_size*60*4)
    cbuf = StringIO()
    cbuf.write('\n'.join([f"{ra} {dec}"]))
    cbuf.seek(0)
    ps1filename = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    r = requests.post(ps1filename, data=dict(filters=filters, type="stack"),files=dict(file=cbuf))
    r.raise_for_status()
    tab = Table.read(r.text, format="ascii")
    urlbase = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    url = urlbase + f"?ra={ra}&dec={dec}&size={size}&format=fits&red={tab['filename'][0]}"
    logger.info(f"Generated Pan-STARRS1 url as {url}")
    return url

def PS1_retrieve(pos, side):
    """get and process image from PS1 using generated url from PS1_geturl"""
    coords = astropy.coordinates.get_icrs_coordinates(pos).to_string(sep=':').split(' ')
    ra, dec = float(coords[0]), float(coords[1])
    fitsurl = PS1_geturl(ra, dec, side, filters='i')
    urllib.request.urlretrieve(fitsurl, "ps1data2.fits")
    image_file = get_pkg_data_filename('ps1data2.fits')
    img = fits.getdata(image_file, ext=0)
    fit = fits.open(fitsurl)[0]
    header = fit.header
    img[np.isnan(img)] = 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wcsf = wcs.WCS(header)
    transform = AsinhStretch() + PercentileInterval(99.75)
    img = transform(img)
    os.remove(str(os.getcwd()) + "/ps1data2.fits")
    logger.info(f"Retrieved image from Pan-STARRS1.")
    return img, header, wcsf

def get_info(gaia_results, center):
    """retrieve additional information from Gaia and Simbad for info box"""
    try:
        gaiaPos = gaia_results['ra'][0], gaia_results['dec'][0]
    except IndexError:
        gaiaPos = center
        logger.info(f"No Gaia sources found in image region.")
        return "NoID", gaiaPos, "None", "None", "None"
    gaia_objID = gaia_results['source_id'][0]
    location = (round(gaiaPos[0], 5), round(gaiaPos[1], 5))
    loc = "Center of search is " + str(location)
    plt.suptitle(f"Gaia ID: {gaia_objID}", fontsize=30, color='black', y=0.95)
    cSimbad = Simbad()
    bibcodes = 'bibcodelist(2000-2030)'
    cSimbad.add_votable_fields('sp', 'pm', bibcodes, 'typed_id')
    cSimbad.ROW_LIMIT = 1 
    ra_r, dec_d = str(gaiaPos[0]), str(gaiaPos[1])
    position = SkyCoord(ra_r, dec_d, unit='deg', frame='icrs')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = cSimbad.query_region(position, radius='0.1m')
    if not results:
        logger.warning(f"No Simbad data found at target ({ra_r, dec_d}).")
        return gaia_objID, gaiaPos, "noData", "noData", "noData"
    try:
        sptype = results['SP_TYPE'][0]
    except ValueError:
        sptype = 'None'
    try:
        bibcode = results['COO_BIBCODE'][0]
    except ValueError: 
        bibcode = 'None'
    logger.info(f"Simbad data found at {ra_r, dec_d}.")
    return gaia_objID, gaiaPos, results, sptype, bibcode

def info_box(gaiaPos, results, sptype, bibcode, survey, side):
    """add info box to bottom left of image"""
    try:
        info = f"Postion: ({round(gaiaPos[0], 4)}, {round(gaiaPos[1], 4)}) \nSurvey: {survey} \nType: {sptype} \nPM: {results['PMRA'][0], results['PMDEC'][0]} \nRef: {bibcode}"
    except TypeError:
        try:
            info = f"Postion: ({round(gaiaPos[0], 4)}, {round(gaiaPos[1], 4)}) \nSurvey: {survey}"
        except TypeError:
            info = f"Position: ({center}) \nSurvey: {survey}"
    props = dict(boxstyle='square', facecolor='white', alpha=0.8)
    ax.annotate(info, xy=(0, 0), xytext=(10, 10), fontsize=10, xycoords='axes fraction', textcoords='offset points', bbox=props, horizontalalignment='left', verticalalignment='bottom')

def plot_setup(wcsf):
    """sets up image plot"""
    global fig
    global ax
    fig = plt.figure(figsize=(15, 15))
    ax = plt.subplot(projection=wcsf)
    fig.tight_layout
    ax.set_autoscale_on(True)
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='white', ls='dotted', linewidth=1)
    
def plot_show_save(positions, plt, fig, img, fname):
    plt.imshow(img, cmap='gray')
    plt.savefig(fname)
    fig.clear()

def main(positions, surveys, side, markers=True, useAllSurveys=False, fileEntry=0, redundant=False):
    if fileEntry!=0:
        positions = fileEntryParse(positions, fileEntry)
    logger.warning(f"Creating findercharts for {len(positions)} position(s).")
    for count, pos in enumerate(positions):
        if fileEntry==0:
            pos = pos.replace(" ", ", ")            
        for num, survey in enumerate(surveys):
            logger.info(f"Creating finderchart {count+1} of {len(positions)} at ({pos}) using {survey}...")
            gaia_objID = ""
            try:
                if "STARR".casefold() in survey.casefold() or "ps1".casefold() in survey.casefold():
                    survey = "Pan-STARRS1"
                    img, header, wcsf = PS1_retrieve(pos, side)
                else:
                    pos2, img, header, wcsf = getImage(pos, survey, side)
            except ConnectionError as e: 
                logger.error("No image found; will be skipping to next image/survey (if applicable).")
                if redundant==False:
                    break
            except astropy.coordinates.name_resolve.NameResolveError as er: 
                if fileEntry==0:
                    logger.error(f"{pos} couldn't be resolved by astropy as a coordinate. If using file entry, be sure to include -f/--file argument.")
                else:
                    logger.error(f"{pos} couldn't be resolved by astropy as a coordinate. Skipping to next image/survey (if applicable).")
            else:
                plot_setup(wcsf)
                center, region = makeGaiaRegion(pos, img, wcsf, side)
                try:
                    gaia_results, gaia_df = regionGaiaQuery(center, region)
                except ConnectionError as e:
                    plot_show_save(positions, plt, fig, img, os.getcwd() + f"/{survey}_noGaiaData.png")
                    logger.error("No Gaia data found for region.")
                    logger.info(f"Finderchart {count+1} created and saved with no Gaia data.\n")
                else: 
                    if markers==True:
                        trans_results = epochAdjust(gaia_results, header, survey)
                        plotPointandPM(trans_results)
                        plotBinaries_fast(gaia_df)
                        gaia_objID, gaiaPos, results, sptype, bibcode = get_info(gaia_results, center)
                        info_box(gaiaPos, results, sptype, bibcode, survey, side)
                    if markers==False:
                        gaia_objID, gaiaPos, results, sptype, bibcode = get_info(gaia_results, center)
                        gaia_objID = str(gaia_objID)
                        gaia_objID += "_noMarks"
                    plot_show_save(positions, plt, fig, img, os.getcwd() + f"/{survey}_{gaia_objID}.png")
                    logger.warning(f"Finderchart {count+1} created and saved.\n")
                    if redundant==False:
                        break
            finally:
                if useAllSurveys==False and redundant==False:
                    break

##/////////////////INPUTS////////////////////////////////////////////
##position resolvable by SIMBAD
##input as list of str
##include 'help' or 'Help' in survey for a list of available surveys. "ps1" or "panstarr" will search PANSTARRS, non-case sensitive.
##input as list of str
#pos = ['1.895920851, -16.09240067']
pos=['300.77531563789, -8.40670962104']
pos = ['118.544, 12.6626']
#pos='Table_spec_params_copy.csv'
positions = ['m83', 'm99']
survey = ['SDSSg', 'ps1']
##image side length in arcmin
side = 0.5
markers = True
useAllSurveys = False
fileEntry = 0
redundant = True
##////////////////END INPUTS/////////////////////////////////////////
if __name__ == '__main__':
    start=time.perf_counter()
    logger, stdlog = log_setup()
    try:
        args = parser()
    except SystemExit as e:
        if '2' in str(e):
            logger.info("Running from Python editor.")
        else:
            raise e
    else:
        if args.quiet==True:
            logger.setLevel(logging.WARNING)
            stdlog.setLevel(logging.WARNING)
        positions = ' '.join(args.pos)
        positions = positions.replace(', ', ',').split(',')
        logger.info(f"Processed positions are {positions}.")
        survey = [*args.survey]
        side, markers, useAllSurveys, fileEntry, redundant = args.width, args.no_markers, args.use_all_surveys, args.file_entry, args.survey_redundant
        if fileEntry!=0:
            positions = args.pos[0]
    main(positions, survey, side, markers, useAllSurveys, fileEntry, redundant)
    end=time.perf_counter()
    logger.info(f"Runtime was {round(end-start, 3)} seconds.")


# In[ ]:




