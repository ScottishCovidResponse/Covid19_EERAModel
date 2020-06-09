import os
import re

class HTMLFileBuilder(object):
    def __init__(self, name="SCRC: COVID-19 EERA Model"):
        self._data = {}
        self._name = name
        self._input_file = None
        self._pages = {'cppcheck' : {'source' : 'cppcheck/index.html',
                             'title' : 'CPP Check',
                             'filename' : 'cppcheck.html'},
                      'doxygen' : {'source' : '../doxygen/html/index.html',
                                  'title' : 'Doxygen Code Documentation',
                                  'filename' : 'doxygen-docs.html'},
                      'coverage' : {'source' : '../build/Covid19EERAModel_coverage/index.html',
                                    'title' : 'Code Coverage',
                                    'filename' : 'code-coverage.html'},
                      'documentation' : {'source' : '../doc/eera_model_overview.html',
                                          'title' : 'Model Overview',
                                          'filename' : 'model_documentation.html'}}

    def build_index(self):
      _index_str='''
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
          <a class="navbar-brand" href="index.html">COVID-19 EERA Model</a>
        
          <!-- Links -->
          <ul class="navbar-nav">
            <!-- Dropdown -->
            <li class="nav-item dropdown">
              <a class="nav-link dropdown-toggle" href="#" id="navbardrop" data-toggle="dropdown">
                Documentation
              </a>
              <div class="dropdown-menu">
                <a class="dropdown-item" href="site/model_documentation.html">Model Overview</a>
                <a class="dropdown-item" href="site/doxygen-docs.html">Doxygen</a>
              </div>
            </li>

            <!-- Dropdown -->
            <li class="nav-item dropdown">
              <a class="nav-link dropdown-toggle" href="#" id="navbardrop" data-toggle="dropdown">
                Code Checks
              </a>
              <div class="dropdown-menu">
                <a class="dropdown-item" href="site/code-coverage.html">Code Coverage</a>
                <a class="dropdown-item" href="site/clang_tidy.html">Clang Tidy</a>
                <a class="dropdown-item" href="site/cppcheck.html">CPP Check</a>
              </div>
            </li>
          </ul>
        </nav> 
            
    </head>
    <body>
        <div class="container">
            <div class="jumbotron">
            <h1>COVID-19 EERA Model</h1>      
            <p>This is the documentation website for the EERA COVID EERA Model repository on GitHub.</p>
            </div>
        </div>
        <div class="container">
          <h2>Documentation</h2>
            <a href="site/model_documentation.html" class="btn btn-info" role="button">Model Overview</a>
            <a href="site/doxygen-docs.html" class="btn btn-info" role="button">Doxygen Documentation Site</a>     
         </div>
        <div class="container">
          <h2>Code Check Reports</h2>   
            <a href="site/code-coverage.html" class="btn btn-primary" role="button">Code Coverage Report</a>   
            <a href="site/clang_tidy.html" class="btn btn-primary" role="button">Clang Tidy Report</a>
            <a href="site/cppcheck.html" class="btn btn-primary" role="button">CPP Check</a>
         </div>
        <div class="container">
            <p></p>
            <h2>Builds</h2>           
            <table class="table table-striped">
              <thead>
                <tr>
                  <th>Branch</th>
                  <th>Status</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td>master</td>
                  <td> <img src="https://github.com/ScottishCovidResponse/Covid19_EERAModel/workflows/Covid19EERAModel/badge.svg?branch=master&event=push" class="img-rounded" alt="Master Status"> </td>
                </tr>
                <tr>
                    <td>dev</td>
                    <td> <img src="https://github.com/ScottishCovidResponse/Covid19_EERAModel/workflows/Covid19EERAModel/badge.svg?branch=dev&event=push" class="img-rounded" alt="Dev Status"> </td>
                </tr>
                <tr>
                    <td>gh-pages</td>
                    <td> <img src="https://github.com/ScottishCovidResponse/Covid19_EERAModel/workflows/Covid19EERAModel/badge.svg?branch=gh-pages&event=push" class="img-rounded" alt="GitHub Pages Status"> </td>
                </tr>
              </tbody>
            </table>
          </div>
          
    </body>
</html>
      '''.format(site_name=self._name)

      with open('index.html', 'w') as f:
        f.write(_index_str)

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

    def build_wrapper_pages(self):
      
      template='''
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
          <a class="dropdown-item" href="{docs_addr}">Model Overview</a>
          <a class="dropdown-item" href="{doxy_addr}">Doxygen</a>
        </div>
      </li>
        <!-- Dropdown -->
        <li class="nav-item dropdown">
          <a class="nav-link dropdown-toggle" href="../index.html" id="navbardrop" data-toggle="dropdown">
            Code Checks
          </a>
          <div class="dropdown-menu">
            <a class="dropdown-item" href="{cov_addr}">Code Coverage</a>
            <a class="dropdown-item" href="{clang_addr}">Clang Tidy</a>
            <a class="dropdown-item" href="{cppc_addr}">CPP Check</a>
          </div>
        </li>
      </ul>
    </nav> 
        
</head>
    <body>
        <div class="container">
            <div class="jumbotron">
            <h1>{title}</h1>      
            <p></p>
            </div>
            <div>
                <center><object data="{address_of_source}" width=1100 height=900></object></center>
            </div>
        </div>
    </body>
</html>
      '''

      for page in self._pages:
          with open('site/'+self._pages[page]['filename'], 'w') as f:
            f.write(template.format(title=self._pages[page]['title'], 
                                    docs_addr=os.path.join('..', 'site', self._pages['documentation']['filename']),
                                    doxy_addr=os.path.join('..', 'site', self._pages['doxygen']['filename']),
                                    cov_addr=os.path.join('..', 'site',self._pages['coverage']['filename']),
                                    clang_addr=os.path.join('..', 'site', self._pages['cppcheck']['filename']).replace('cppcheck','clang_tidy'),
                                    cppc_addr=os.path.join('..', 'site', self._pages['cppcheck']['filename']),
                                    address_of_source=self._pages[page]['source'], site_name=self._name))

    def build_clang_html(self):
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
        
        with open('site/clang_tidy.html', 'w') as f:
          f.write(_body.format(body=_table_body, site_name=self._name))
    
    def Run(self):
      self.build_index()
      self.build_wrapper_pages()
      self.build_clang_html()
        

if __name__ in "__main__":
    import argparse

    parser = argparse.ArgumentParser('ClangTidyParser')

    parser.add_argument('input_file', help='Build log file')

    cl_parser = HTMLFileBuilder()
    cl_parser.parse(parser.parse_args().input_file)
    cl_parser.Run()
    
    
