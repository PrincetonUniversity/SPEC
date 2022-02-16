window.MathJax = {
  loader: {
    load: ['[tex]/color', '[tex]/cancel']
  },
  tex: {
      tags: 'ams',
      packages: {'[+]': ['noerrors'],
                 '[+]': ['color'],
                 '[+]': ['verb'],
                 '[+]': ['cancel']
      },
      macros: { bm: ['{\\boldsymbol #1}',1]
      }
  }
};

// https://docs.mathjax.org/en/latest/web/configuration.html#configuring-and-loading-in-one-script
//(function () {
//  var script = document.createElement('script');
//  script.src = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js';
//  script.async = true;
//  document.head.appendChild(script);
//})();
