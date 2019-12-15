MathJax.Hub.Config({
    TeX: { // Enable automatic equation numbering
        equationNumbers: { autoNumber: "AMS" }
    }
});
MathJax.Hub.Config({ TeX: { extensions: ["color.js"] }});
MathJax.Hub.Register.StartupHook("TeX color Ready", function() {
    // https://stackoverflow.com/questions/44715432/how-to-define-custom-colors-in-mathjax-in-config
    var color = MathJax.Extension["TeX/color"];
    color.colors["red"] = color.getColor('RGB','255,0,0');
    color.colors["gre"] = color.getColor('RGB','0,255,0');
    color.colors["blu"] = color.getColor('RGB','0,0,255');
    color.colors["ora"] = color.getColor('RGB','255,127,0');
    color.colors["gra"] = color.getColor('RGB','0,127,255');
});
