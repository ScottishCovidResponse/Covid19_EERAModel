
import re

class ClangTidyParser(object):
    def __init__(self):
        self._data = {}
        self._input_file = None

    def parse(self, input_file):
        self._input_file = input_file
        with open(self._input_file) as f:
            for line in f:
                try:
                    entry = re.findall('.+(src/\w+[.h|.cpp].+)', line)[0]
                    entry  = entry.split(': ', 1)
                    try:
                        file_name, address = entry[0].split(':', 1)
                    except ValueError:
                        file_name = entry[0]
                        address = ''
                    library = re.findall('\[(.+)\]', entry[1])
                    message = entry[1] if len(library) == 0 else entry[1].replace(' [{}]'.format(library[0]), '')
                    if file_name not in self._data:
                        self._data[file_name] = {}
                    self._data[file_name][address] = {'message' : message, 
                                                 'library' : '' if len(library) == 0 else library[0]}
                except IndexError:
                    continue
    
    def build_html(self, site_name="SCRC: COVID-19 EERA Model"):
        _body='''
<!DOCTYPE HTML>
<html lang="en-GB">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.5.0/css/bootstrap.min.css">
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js"></script>
        
        <title>{site_name}</title>
         <!-- A grey horizontal navbar that becomes vertical on small screens -->
        <nav class="navbar navbar-expand-sm bg-dark navbar-dark">

        <!-- Links -->
        <ul class="navbar-nav">
            <li class="nav-item">
            <a class="nav-link" href="../index.html">Home</a>
            </li>
            <li class="nav-item">
            <a class="nav-link" href="../doxygen/html/index.html">Doxygen</a>
            </li>
            <li class="nav-item">
            <a class="nav-link" href="../build/Covid19EERAModel_coverage/index.html">Code Coverage</a>
            </li>
            <li class="nav-item">
            <a class="nav-link" href="doxygen/clang_tidy.html">Clang Tidy</a>
            </li>
        </ul>

        </nav> 
            
    </head>
    <body>
    <div class="container">
    <div class="jumbotron">
        <h1>Clang Tidy</h1>
        <p>Listed below are the warnings and suggestions made by Clang Tidy during compilation.</p>
    </div>
    <p></p>
    </div>
    <div class="container">
    <table class="table table-striped">
    <thead>
      <tr>
        <th>Type</th>
        <th>Address</th>
        <th>Library</th>
        <th>Message</th>
      </tr>
    </thead>
    {body}
    </container>
    </body>
</html>
'''
        _buttons = {'warning' : '<button type="button" class="btn btn-warning"><strong>Warning</strong></button>',
                    'info' : '<button type="button" class="btn btn-info"><strong>Suggestion</strong></button>'}

        _table_body = ''

        for entry in self._data:
            entry_str='''
            <tr class="table-dark text-dark">
            <td colspan="4" style="text-align:center"><strong>{file}</strong></td>
            </tr>
            '''.format(file=entry)
            for address in sorted(self._data[entry].keys()):
                is_warning = 'warning' in self._data[entry][address]['message'].lower()
                button = _buttons['warning'] if is_warning else _buttons['info']

                entry_str+='''
                <tr>
                    <td>{button}</td>
                    <td>{address}</td>
                    <td>{library}</td>
                    <td>{message}</td>
                </tr>
                '''.format(button=button, file=entry, address=address,
                            library=self._data[entry][address]['library'], message=self._data[entry][address]['message'])
            _table_body += entry_str
        
        return _body.format(body=_table_body, site_name=site_name)
        

if __name__ in "__main__":
    import argparse

    parser = argparse.ArgumentParser('ClangTidyParser')

    parser.add_argument('input_file', help='Build log file')

    cl_parser = ClangTidyParser()
    cl_parser.parse(parser.parse_args().input_file)
    
    with open('doxygen/clang_tidy.html', 'w') as f:
        f.write(cl_parser.build_html())
