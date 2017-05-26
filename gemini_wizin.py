#!/usr/bin/env python
import os
import warnings
import webbrowser
import subprocess


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
                        " {skip_gerp_bp} {skip_cadd}"
                        "{cores} {dbname}.db")
    return gemini_load_cmd




@app.route('/', methods='POST')
def wizin():
    # bottom
    load = request.GET.get('load')

    # files and options
    vcf = str(request.GET.get('VCFfile','').strip())
    dbname = str(request.GET.get('outfilename','').strip())

    cnv = request.GET.get('cnv','').strip()
    if cnv =='on':
        cnv = '--cnv'

    CNVmap = str(request.GET.get('CNVmap','').strip())
    if CNVmap:
        CNVmap = '--dgv_cnvmap ' + CNVmap

    ped_file = str(request.GET.get('PED','').strip())
    if ped_file:
        ped_file = '-p ' + ped_file

    skip_gerp_bp = str(request.GET.get('skip_gerp_bp','').strip())
    if skip_gerp_bp =='on': skip_gerp_bp = '--skip-gerp-bp'

    skip_cadd = str(request.GET.get('skip_cadd','').strip())
    if skip_cadd =='on': skip_cadd = '--skip-cadd'

    cores = str(request.GET.get('CORES','').strip())
    if cores:
        cores = '--cores ' + cores

    submit_command = "{cmd}"

    if load=="True":
        gemini_load_c = gemini_load_wiz().format(**locals())
        print gemini_load_c
        subprocess.call(gemini_load_c, shell = True)
        return template('wizin.j2')
    else:
        return template('wizin.j2')



def wizin_main(parser, args):
    host = args.host
    port = args.port

    webbrowser.open_new_tab("http://{}:{}".format(host, port))
    run(app, host=host, port=port,
        reloader=True, debug=True)

    try:
        if args.use == "puzzle":
            browser_puzzle(args)
            # XXX: https://github.com/dgaston/kvasir
            #if args.use == "kvasir":
            #    raise NotImplementedError
        elif args.use == "builtin":
            browser_builtin(args)
        else:
            raise NotImplementedError("GEMINI-compatible Browser '{browser}' not found.".format(browser=browser))
    except ImportError as e:
        raise ImportError("{exc}\nIs {browser} correctly installed?".format(exc=e, browser=browser))
