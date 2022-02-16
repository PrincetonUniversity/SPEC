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
