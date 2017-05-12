import os
import webbrowser
import subprocess

# -- common bottle importation
from bottle import TEMPLATE_PATH, Bottle, run, static_file, debug, request
from bottle import jinja2_template as template
debug(True)

base_dir = os.path.dirname(__file__)
TEMPLATE_PATH.append(os.path.abspath(os.path.join(base_dir, 'views')))

# -- the instance app is important
app = Bottle()

# -- serve static files, files located in static
static_folder = 'static'
_static_folder = os.path.join(os.path.dirname(__file__), static_folder)

@app.route('/wiz', method='POST')
@app.route('/')
def wizin():
    if request.GET.get('load', '').strip():
        VCFfile = str(request.GET.get('CNVmap').strip())
        #cnv = request.GET.get('cnv')
        #CNVmap = str(request.GET.get('CNVmap').strip())
        outfile = str(request.GET.get('outfilename').strip())

        gemini_load_cmd = ("gemini_cnv load -v " + VCFfile + " " + outfile +".db")

        subprocess.check_call(gemini_load_cmd, shell=True, stderr=sys.stderr)
    else:
        return template('wizin.j2')


def wizin_main(parser, args):
    host = args.host
    port = args.port

    webbrowser.open_new_tab("http://{}:{}".format(host, port))
    run(app, host=host, port=port,
        reloader=True, debug=True)
