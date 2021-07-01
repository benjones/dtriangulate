module dtriangulate.svg;

import std.stdio : File;

//write the header and return points scaled appropriately
Vec[] prepareSVG(Vec)(ref File f, const Vec[] points, const size_t[] activePoints){
  import std.array;
  import std.stdio : File;
  import std.algorithm;

  //activePoints.sort().uniq();
  if(activePoints.empty){
	f.writeln("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1200\" height=\"1200\" ></svg>");
	return [];
  }
  //using a delecate to store a copy of the points, not a reference
  Vec[] usedPoints = activePoints.map!(delegate Vec(size_t a){ return points[a]; }).array;


  Vec minP = Vec(usedPoints.map!("a.x").minElement(), usedPoints.map!("a.y").minElement());
  Vec maxP = Vec(usedPoints.map!("a.x").maxElement(), usedPoints.map!("a.y").maxElement());

  Vec center = 0.5*(minP + maxP);

  Vec size = maxP - minP;
  minP -= .03*size;
  maxP += .03*size;

  size = maxP - minP;
  auto maxDim = max(size.x, size.y);

  auto radius = max(size.x, size.y)/400.0f;

  auto aspectRatio = size.y/size.x;

  auto width = 2000;
  auto height = 2000;//*aspectRatio;


  Vec[] scaledPoints = points.map!( a => (a - center)*(20/maxDim))
	//	.map!(a => Vec(a.x*2.5, a.y))
	.array;

  f.write("<svg xmlns=\"http://www.w3.org/2000/svg\"  ");
  f.writefln("width=\"%d\" height=\"%d\" viewBox=\"-10 -10 20 20\" >", width, height);

  f.writeln("<g transform=\"scale(1, -1)\" >");

  return scaledPoints;

}

void svgLine(Vec)(ref File f, Vec a, Vec b, float lineWidth, bool dashed = false, string color = "black"){
  string dashedString = dashed ? "stroke-dasharray=\"2%, 10%\"" : "";
  f.writefln("<line x1=\"%.8f\" y1=\"%.8f\" x2=\"%.8f\" y2=\"%.8f\" stroke-width=\"%.8f\" %s stroke=\"%s\"  />",
             a.x, a.y, b.x, b.y, lineWidth, dashedString, color);

}

void svgLabeledPoints(Vec)(ref File f, const Vec[] scaledPoints, const size_t[] activePoints){

  foreach(i ; activePoints){
	const auto p = scaledPoints[i];
	f.writefln( "<circle cx=\"%.8f\" cy=\"%.8f\" r=\"%.8f\" fill=\"black\" />", p.x, p.y, .1);
	f.writefln( "<g transform=\"translate(%.8f, %.8f) scale(1, -1)\" >" ~
				"<text x=\"0\" y=\"0\" font-family=\"Verdana\" font-size=\"%.8f\" fill=\"blue\" >%d</text></g>",
				p.x, p.y, .45, i);
  }
}
