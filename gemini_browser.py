import os
import warnings
import webbrowser
import subprocess
import shutil
from collections import namedtuple
import numpy as np

import GeminiQuery

import GeminiQuery
from gemini.gim import (AutoDom, AutoRec, DeNovo, MendelViolations, CompoundHet)

database = None

# based upon bottle example here:
# https://bitbucket.org/timtan/bottlepy-in-real-case

# -- determine where I launch python and config lib path
# base_dir = os.path.dirname(__file__)
# third_party_path = os.path.abspath(os.path.join(base_dir, 'third_party' ))
# sys.path.insert(0, third_party_path)

# -- common bottle importation
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from bottle import TEMPLATE_PATH, Bottle, run, static_file, debug, request, redirect, FileUpload
    from bottle import jinja2_template as template

debug(True)


base_dir = os.path.dirname(__file__)
TEMPLATE_PATH.append(os.path.abspath(os.path.join(base_dir, 'views')))

# -- the instance app is important
app = Bottle()

# -- serve static files, files located in static
static_folder = 'static'
_static_folder = os.path.join(os.path.dirname(__file__), static_folder)

def name_sample(database):
	samples = []
	rows_sample_query = 'select name, paternal_id, maternal_id from samples'
	rows_sample = GeminiQuery.GeminiQuery(database)
	rows_sample._set_gemini_browser(True)
	rows_sample.run(rows_sample_query)
	for row in rows_sample:
		samples.append(row)
	return samples

@app.route('/stats/region/:chrom', method='GET')
def stats_region(chrom):
    # Note: chrom is give as an argument

    # we then extract start and end using HTML GET
    start = request.GET.get('start', '').strip()
    end = request.GET.get('end', '').strip()

    # construct a query
    query =  "SELECT start, end from variants"
    query += " WHERE chrom = '" + chrom + "'"
    query += " AND start >= " + start
    query += " AND end <= " + end

    # issue the query
    gq = GeminiQuery.GeminiQuery(database)
    gq._set_gemini_browser(True)
    gq.run(query)

    # return query results in JSON format
    return{'features': [dict(row) for row in gq]}


@app.route('/static/<filepath:path>')
def server_static(filepath):
    return static_file(filepath, root=_static_folder)
# -- end of static folder configuration

# -- index page routing


### loading wizard ###
@app.route('/index')
@app.route('/')
def index():
    # control the exist of database in the path make sure database is found if provided
    if database is None or not os.path.exists(database):
        redirect('/wizin')

    return template('index.j2')

@app.route('/query_json', method='GET')
def query_json():
    query = request.GET.get('query', '').strip()

    gq = GeminiQuery.GeminiQuery(database)
    gq._set_gemini_browser(True)
    gq.run(query)

    return {'gemini_results': [dict(row) for row in gq]}


@app.route('/query', method='GET')
def query():

    def _get_fields():
        query = request.GET.get('query', '').strip()
        gt_filter = request.GET.get('gt_filter', '').strip()
        use_header = request.GET.get('use_header')
        igv_links = request.GET.get('igv_links')
        return query, gt_filter, use_header, igv_links

    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():

        (query, gt_filter, use_header, igv_links) = _get_fields()

        if use_header: use_header = True
        if igv_links: igv_links = True

        gq = GeminiQuery.GeminiQuery(database)
        gq._set_gemini_browser(True)
        gq.run(query, gt_filter)


        if len(query) == 0:
            return template('query.j2', dbfile=database)

        if igv_links and ('chrom' not in query.lower()
                          or 'start' not in query.lower()
                          or 'end' not in query.lower()):
            return template('query.j2', dbfile=database,
                            rows=gq,
                            igv_links=igv_links,
                            igv_links_error=True,
                            use_header=use_header,
                            gt_filter=gt_filter,
                            query=query)
        else:
            return template('query.j2', dbfile=database,
                            rows=gq,
                            igv_links=igv_links,
                            igv_links_error=False,
                            use_header=use_header,
                            gt_filter=gt_filter,
                            query=query)

    # user clicked the "save to file" button
    elif request.GET.get('save', '').strip():

        (query, gt_filter, use_header, igv_links) = _get_fields()

        gq = GeminiQuery.GeminiQuery(database)
        gq.run(query, gt_filter)

        if len(query) == 0:
            return template('query.j2', dbfile=database)

        # dump the results to a text file.  this will be
        # stored in /static and a link will be given to
        # the user.
        tmp_file = '/tmp.txt'
        tmp = open(_static_folder + tmp_file, 'w')
        for row in gq:
            tmp.write('\t'.join(str(c) for c in row) + '\n')
        tmp.close()

        return template('query.j2', dbfile=database,
                        tmp_file=tmp_file,
                        igv_links=igv_links,
                        igv_links_error=True,
                        use_header=use_header,
                        gt_filter=gt_filter,
                        query=query)
    # user did nothing.
    else:
        return template('query.j2', dbfile=database)


default_cols = ['chrom', 'start', 'end', 'ref', 'alt',
                'polyphen_pred', 'sift_pred',
                'max_aaf_all', 'impact', 'impact_severity',
                'gene', 'biotype', 'transcript']
# turn a dictionary into something that can be accessed by attribute.
class Arguments(object):
    """
    >>> args = Arguments(db='some.db')
    >>> args.db
    'some.db'
    """
    _opts = ("columns", "db", "filter", "min_kindreds", "families",
                 "pattern_only", "max_priority", # only for compound_het
                 "allow_unaffected", "min_gq", "lenient", "min_sample_depth"
                 "f" # only for overlap
                 )
    def __init__(self, **kwargs):
        if not "min_gq" in kwargs: kwargs['min_gq'] = 0
        if not "lenient" in kwargs: kwargs['lenient'] = False
        for k in ("families", "filter"):
            if not k in kwargs: kwargs[k] = None
        if not "gt_phred_ll" in kwargs: kwargs['gt_phred_ll'] = None
        if not "min_sample_depth" in kwargs: kwargs['min_sample_depth'] = 0

        for k in ("min_kindreds", "max_priority"):
            if not k in kwargs: kwargs[k] = 1
        for k in ("pattern_only", "allow_unaffected"):
            if not k in kwargs: kwargs[k] = False
        if not "columns" in kwargs:
            kwargs['columns'] = ",".join(default_cols)
        if not "f_par" in kwargs: kwargs['f_par'] = None
        if not 'r' in kwargs: kwargs['r']=None
        if not 'alt_par' in kwargs: kwargs['alt_par'] = None
        if not 'int_len_min' in kwargs: kwargs['int_len_min'] = None
        if not 'int_len_max' in kwargs: kwargs['int_len_max'] = None
        if not 'gene_map' in kwargs: kwargs['gene_map'] = None
        if not 'sample' in kwargs: kwargs['sample'] = None
        if not 'v' in kwargs: kwargs['v'] = None
        if not 'dgv_cnvmap' in kwargs: kwargs['dgv_cnvmap'] = None
        self.__dict__.update(**kwargs)



@app.route('/de_novo', method='GET')
def de_novo():

    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():

        min_sample_depth = str(request.GET.get('min-depth', '10').strip())
        igv_links = request.GET.get('igv_links')

        args = Arguments(db=database, min_sample_depth=min_sample_depth)

        row_iter = DeNovo(args).report_candidates()

        return template('de_novo.j2', dbfile=database,
                        rows=row_iter,
                        igv_links=igv_links)

    else:
        return template('de_novo.j2', dbfile=database)


@app.route('/auto_rec', method='GET')
def auto_rec():

    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():
        min_sample_depth = str(request.GET.get('min-depth', '10').strip())
        args = Arguments(db=database, min_sample_depth=min_sample_depth)

        row_iter = AutoRec(args).report_candidates()
        return template('auto_rec.j2', dbfile=database, rows=row_iter)

    else:
        return template('auto_rec.j2', dbfile=database)


@app.route('/auto_dom', method='GET')
def auto_dom():

    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():
        min_sample_depth = str(request.GET.get('min-depth', '10').strip())
        args = Arguments(db=database, min_sample_depth=min_sample_depth)
        row_iter = AutoDom(args).report_candidates()
        return template('auto_dom.j2', dbfile=database, rows=row_iter)

    else:
        return template('auto_dom.j2', dbfile=database)


@app.route('/db_schema', method='GET')
def db_schema():
    return template('db_schema.j2')


@app.route('/overlap', method=['POST','GET'])
def overlap():
    import tool_overlap

    name = tool_overlap.name_dgv(database=database)
    if name == None:
        name = 'no map'

    rows_sample = name_sample(database=database)

    f_str = request.POST.get('f')
    recip = request.POST.get('reciprocal')
    len_min = request.POST.get('int_len_min')
    len_max = request.POST.get('int_len_max')
    alt = request.POST.get('alt_par')
    sample = request.POST.get('sample')
    invert = request.POST.get('invert')
    dgv_cnvmap = request.files.get('CNVmap')

    args = Arguments(db=database, f_par = f_str, int_len_min = len_min, int_len_max = len_max, alt_par=alt,sample = sample, v = invert, dgv_cnvmap = dgv_cnvmap)


    if dgv_cnvmap:
        name, ext = os.path.splitext(dgv_cnvmap.filename)
        # control extension
        if ext !='.txt':
            message = 'File extension of CVN map is not allowed.'
            return template('overlap.j2', message=message)

        # save file
        save_path = os.getcwd()
        file_path = "{path}/{file}".format(path=save_path, file = dgv_cnvmap.raw_filename)
        dgv_cnvmap.save(file_path, overwrite = True)
        args.dgv_cnvmap = dgv_cnvmap.raw_filename
        name = dgv_cnvmap.raw_filename

    if recip != None:
        args.r = True
    if invert != None:
        args.v = True

    # user clicked the "submit" button
    if request.POST.get('submit', '').strip():
        tool_overlap.overlap_res(args)
        if args.v == True:
            query_all = "SELECT * FROM no_overlap"
        else:
            query_all = "SELECT * FROM overlap"
        over = GeminiQuery.GeminiQuery(args.db)
        over._set_gemini_browser(True)
        over.run(query_all)

        result = 'Results with: '
        if f_str != '': result += 'minimum overlap fraction = ' + f_str
        if recip != None: result += ' , reciprocal = ' + recip
        if invert != None: result += ' , no overlap result = ' + invert
        if len_min != '': result += ' , minimum overlap length = ' + len_min
        if len_max != '': result += ' , maximum overlap length = ' + len_max
        if alt != '': result += ', alteration = ' + alt
        else: result += ' -'
        return template('overlap.j2', dbfile=database, rows=over, maps_name = name, results = result, reciprocal = recip, rows_sample = rows_sample, invert = invert)

    # user clicked the "save as a text file" button
    elif request.POST.get('save', '').strip():
        res = request.POST.get('results')

        tmp_file = 'overlap_result.txt'
        tmp = open(tmp_file, 'w')
        tmp.write('## ' + res + '\n')

        query_all = "SELECT * from overlap"
        result = GeminiQuery.GeminiQuery(database)
        result._set_gemini_browser(True)
        result.run(query_all)
        tmp.write('#'+ str(result.header) + '\n')
        for row in result:
            tmp.write(str(row)+'\n')
        tmp.close()
        return template('overlap.j2', dbfile=database, maps_name = name, rows_sample = rows_sample, rows=result)
    else:
        return template('overlap.j2', dbfile=database, maps_name = name, rows_sample = rows_sample )

@app.route('/over_gene', method=['POST','GET'])
def overlap_gene():
    import tool_overlap_gene
    alt = np.array([])

    gene_map = request.files.get('genemap')
    gen_check = request.POST.get('gen_check')
    sample = request.POST.get('sample')

    rows_sample = name_sample(database=database)

    args = Arguments(db = database,gene_map = gene_map,sample = sample)
    name_map = 'Ensembl 75 version'

    if gene_map:
        name, ext = os.path.splitext(gene_map.filename)

        # control extension
        if ext not in ('.txt','.bed'):
            message = 'ERR: File extension of gene map file is not allowed.'
            return template('over_gene.j2', message=message)

        #save file
        save_path = os.getcwd()
        file_path = "{path}/{file}".format(path=save_path, file=gene_map.raw_filename)
        gene_map.save(file_path,overwrite=True)
        args.gene_map = gene_map.raw_filename
        name_map = gene_map.raw_filename


    #bottom
    if request.POST.get('submit', '').strip():
        if gene_map != None:
            tool_overlap_gene.get_gene_map(args=args)
            gene,alt,res = tool_overlap_gene.overlap_custom_gene_browser(args)
        else:
            gene, alt,res = tool_overlap_gene.overlap_gene_browser(args)

        return template('over_gene.j2', rows=res, name_map = name_map, gen_check = gen_check, sample = args.sample, rows_sample=rows_sample)


    if request.POST.get('save', '').strip():

        tmp_file = 'overlap_gene_result.txt'
        tmp = open(tmp_file, 'w')

        if gen_check != None:
            gene,alt,res = tool_overlap_gene.overlap_custom_gene_browser(args)
            name_map = "custom map"
        else:
            gene,alt,res = tool_overlap_gene.overlap_gene_browser(args)

        for row in res:
            tmp.write(str(row)+'\n')
        tmp.close()
        return template('over_gene.j2', rows=res, name_map = name_map, gen_check = gen_check, sample = args.sample, rows_sample=rows_sample)

    elif request.POST.get('heatmap','').strip():
        if gen_check != None:
            gene,alt,res = tool_overlap_gene.overlap_custom_gene_browser(args)
        else:
            gene,alt,res = tool_overlap_gene.overlap_gene_browser(args)

        sel_sample = []
        if sample:
            sel_sample = sample.split(',')
        else:
            for r in rows_sample:
                sel_sample.append(r['name'])
        tool_overlap_gene.heatmap(args.db,alt,gene,sel_sample)
        name, ext = str(database).split('.')
    	path_name = os.getcwd() + '/'
        picture = path_name + name + '_gene_heatmap.png'
        webbrowser.open('file://' + picture)
    else:
        return template('over_gene.j2', name_map= name_map, rows_sample = rows_sample)


### loading wizard ###
@app.route('/wizin', method=['POST','GET'])
def wizin():

    def gemini_load_wiz():
        gemini_load_cmd = ("gemini_cnv load -v {vcf} {cnv} {ped_file} "
                            "{cores} %s") %database
        return gemini_load_cmd

    message = "%s does not exists in the current path.<br> Create it with the Loading Wizard below.\n" %database
    # bottom
    load = request.POST.get('load')


    # re-create save folder
    save_path = os.getcwd() + '/' + str(database[:-3])
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    else:
        shutil.rmtree(save_path)
        os.makedirs(save_path)


    ## files and options
    vcf = request.files.get('VCFfile')
    if vcf:
        name, ext = os.path.splitext(vcf.filename)
        # control extension
        if ext not in ('.vcf','.gz'):
            message = 'File extension of VCF file is not allowed.'
            return template('wizin.j2', message=message)

        #save file
        file_path = "{path}/{file}".format(path=save_path, file=vcf.raw_filename)
        vcf.save(file_path)
        vcf = save_path + '/' + vcf.raw_filename


    cnv = request.POST.get('cnv')
    if cnv =='on':cnv = '--cnv'
    else: cnv = ''

    ped_file = request.files.get('PED')
    if ped_file:
        name, ext = os.path.splitext(ped_file.filename)
        # control extension
        if ext !='.ped':
            message = 'File extension of PED file  is not allowed.'
            return template('wizin.j2', message=message)

        # save file
        file_path = "{path}/{file}".format(path=save_path, file=ped_file.raw_filename)
        ped_file.save(file_path)
        ped_file = '-p ' + save_path + '/' + ped_file.raw_filename
    else: ped_file=''

    cores = request.params.get('CORES')
    if cores:
        cores = '--cores ' + cores

    submit_command = "{cmd}"

    if load == 'True':
        gemini_load_c = gemini_load_wiz().format(**locals())
        print gemini_load_c
        subprocess.call(gemini_load_c, shell = True)
        redirect('/index')

    return template('wizin.j2',message=message)



## Switch between the different available browsers
def browser_puzzle(args):
    host = args.host
    port = args.port

    try:
        # Puzzle browser plugin
        from puzzle.server import factory as puzzle_app
        from puzzle.plugins import GeminiPlugin
        from puzzle.server.settings import BaseConfig
    except ImportError as e:
        raise ImportError("%s\nPlease 'pip install puzzle' if you want to run it\n" % e)

    plugin = GeminiPlugin(db=args.db, vtype="sv")
    root = os.path.expanduser("~/.puzzle")

    BaseConfig.PUZZLE_BACKEND = plugin
    BaseConfig.UPLOAD_DIR = os.path.join(root, 'resources')

    puzzle_srv = puzzle_app.create_app(config_obj=BaseConfig)
    webbrowser.open_new_tab("http://{}:{}".format(host, port))
    run(puzzle_srv, host=host, port=port)

def browser_builtin(args):
    host = args.host
    port = args.port

    webbrowser.open_new_tab("http://{}:{}".format(host, port))
    run(app, host=host, port=port,
        reloader=True, debug=True)


def browser_main(parser, args):
    global database
    database = args.db
    browser = args.use

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
