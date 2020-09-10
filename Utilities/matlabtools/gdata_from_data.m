function gdata = gdata_from_data(data)


gdata = data.grid;
gdata.mn = data.output.mn;
gdata.im = data.output.im;
gdata.in = data.output.in;
gdata.Rbc = data.output.Rbc;
gdata.Rbs = data.output.Rbs;
gdata.Zbc = data.output.Zbc;
gdata.Zbs = data.output.Zbs;
gdata.Igeometry = data.input.physics.Igeometry;
gdata.Mregular = data.input.numerics.Mregular;



end