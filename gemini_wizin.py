#!/usr/bin/env python
import os
import sys
import warnings
import webbrowser
import subprocess

# gemini modules
import gemini_browser

# -- common bottle importation
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from bottle import TEMPLATE_PATH, Bottle, run, static_file, debug, request
    from bottle import jinja2_template as template
debug(True)

base_dir = os.path.dirname(__file__)
TEMPLATE_PATH.append(os.path.abspath(os.path.join(base_dir, 'views')))

# -- the instance app
app = Bottle()

# -- serve static files, files located in static
static_folder = 'static'
_static_folder = os.path.join(os.path.dirname(__file__), static_folder)

def gemini_load_wiz():
    gemini_load_cmd = ("gemini_cnv load -v {vcf} {cnv} {CNVmap} {ped_file}"
                        "{cores} {dbname}.db")
    return gemini_load_cmd

def gemini_launch_browser(dbname):
    gemini_launch_browser_cmd = ("gemini_cnv browser " + dbname +".db ")
    return gemini_launch_browser_cmd


@app.route('/', methods='POST')
def wizin():
    # bottom
    load = request.GET.get('load')

    # files and options
    vcf = request.params.get('VCFfile')
    dbname = request.params.get('outfilename')

    cnv = request.GET.get('cnv')
    if cnv =='on':
        cnv = '--cnv'

    CNVmap = request.params.get('CNVmap')
    if CNVmap:
        CNVmap = '--dgv_cnvmap ' + CNVmap

    ped_file = request.params.get('PED')
    if ped_file:
        ped_file = '-p ' + ped_file

    cores = request.params.get('CORES')
    if cores:
        cores = '--cores ' + cores

    submit_command = "{cmd}"

    if load=="True":
        gemini_load_c = gemini_load_wiz().format(**locals())
        print gemini_load_c
        subprocess.call(gemini_load_c, shell = True)
        return template('end.j2')
    else:
        return template('wizin.j2')



def wizin_main(parser, args):
    host = args.host
    port = args.port

    webbrowser.open_new_tab("http://{}:{}".format(host, port))
    run(app, host=host, port=port,
        reloader=True, debug=True)
