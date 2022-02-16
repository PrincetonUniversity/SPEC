// <script type="text/x-mathjax-config">
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

(function () {
  var script = document.createElement('script');
  script.src = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js';
  script.async = true;
  document.head.appendChild(script);
})();

// </script>
// <script type="text/javascript"
//         src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
// </script>
