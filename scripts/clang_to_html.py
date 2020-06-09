
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
    <!-- Redirect -->
    <!--<meta http-equiv="refresh" content="1;url=doxygen/html/index.html">-->
    <meta name="viewport" content="width=device-width, initial-scale=1">
    
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.5.0/css/bootstrap.min.css">
    <script src="https://code.jquery.com/jquery-3.5.1.slim.min.js" integrity="sha384-DfXdz2htPH0lsSSs5nCTpuj/zy4C+OGpamoFVy38MVBnE+IbbVYUew+OrCXaRkfj" crossorigin="anonymous"></script>
    <script src="https://cdn.jsdelivr.net/npm/popper.js@1.16.0/dist/umd/popper.min.js" integrity="sha384-Q6E9RHvbIyZFJoft+2mJbHaEWldlvI9IOYy5n3zV9zzTtmI3UksdQRVvoxMfooAo" crossorigin="anonymous"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/js/bootstrap.min.js" integrity="sha384-OgVRvuATP1z7JjHLkuOU7Xw704+h835Lr+6QL9UvYjZE3Ipu6Tp75j7Bh/kR0JKI" crossorigin="anonymous"></script>
    
    <title>{site_name}</title>

    <nav class="navbar navbar-expand-sm bg-dark navbar-dark">
      <!-- Brand -->
      <a class="navbar-brand" href="../index.html">COVID-19 EERA Model</a>
    
      <!-- Links -->
      <ul class="navbar-nav">
        <!-- Dropdown -->
      <li class="nav-item dropdown">
        <a class="nav-link dropdown-toggle" href="#" id="navbardrop" data-toggle="dropdown">
          Documentation
        </a>
        <div class="dropdown-menu">
          <a class="dropdown-item" href="../site/model_documentation.html">Model Overview</a>
          <a class="dropdown-item" href="../site/doxygen-docs.html">Doxygen</a>
        </div>
      </li>
        <!-- Dropdown -->
        <li class="nav-item dropdown">
          <a class="nav-link dropdown-toggle" href="../index.html" id="navbardrop" data-toggle="dropdown">
            Code Checks
          </a>
          <div class="dropdown-menu">
            <a class="dropdown-item" href="../site/code-coverage.html">Code Coverage</a>
            <a class="dropdown-item" href="../site/clang_tidy.html">Clang Tidy</a>
            <a class="dropdown-item" href="../site/cppcheck.html">CPP Check</a>
          </div>
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
    
    with open('site/clang_tidy.html', 'w') as f:
        f.write(cl_parser.build_html())
