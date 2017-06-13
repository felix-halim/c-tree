(function () {
  var svg = d3.select('#trimmer');

  draw_phase1(svg, 20, 20, 250, 250, 10, 100);
  draw_phase2(svg, 300, 20, 300, 250, 10, 100, 30);
  draw_phase3(svg, 630, 20, 500, 250, 10, 100, 30);

  function draw_phase1(svg, sx, sy, width, height, bwidth, bheightL, bheightS) {
    var g = svg.append('g').attr('transform', 'translate(' + sx + ',' + sy + ')');
    var cx = width / 2, cy = 30;
    // g.append('rect').attr({ width: width, height: height, fill: 'lightgray', });
    g.append('polygon').attr({
      points: cx + ',10 ' + (cx + 10) + ',' + (cy + 5) + ' ' + (cx - 10) + ',' + (cy + 5),
      style: "fill:blue; stroke:black; stroke-width:1"
    });
    g.append('circle').attr({ cx: cx, cy: cy, r: 3, fill: 'blue', });
    draw_chain(g, cx, cy, cx - width / 2 + 20, cy + 50, 7, bwidth, bheightL, 20);
  }

  function draw_phase2(svg, sx, sy, width, height, bwidth, bheightL, bheightS) {
    var g = svg.append('g').attr('transform', 'translate(' + sx + ',' + sy + ')');
    var cx = width / 2, cy = 30;
    // g.append('rect').attr({ width: width, height: height, fill: 'lightgray', });
    g.append('polygon').attr({
      points: cx + ',10 ' + (cx + 50) + ',35 ' + (cx - 50) + ',35',
      style: "fill:blue; stroke:black; stroke-width:1"
    });

    var nextxi = cx - 25, nextx = cx - width / 2 + 20, n = 2, gap = 15;
    draw_chain(g, nextxi, cy, nextx, cy + 50, n, bwidth, bheightL, gap);

    nextxi += 10; nextx += n * (bwidth + gap) + gap; n = 1;
    console.log(nextx);
    draw_chain(g, nextxi, cy, nextx, cy + 50, n, bwidth, bheightL, gap);

    nextxi += 10; nextx += n * (bwidth + gap) + gap; n = 1;
    draw_chain(g, nextxi, cy, nextx, cy + 50, n, bwidth, bheightS, gap);

    nextxi += 10; nextx += n * (bwidth + gap); n = 1;
    draw_chain(g, nextxi, cy, nextx, cy + 50, n, bwidth, bheightS, gap, true);

    nextxi += 10; nextx += n * (bwidth + gap); n = 1;
    draw_chain(g, nextxi, cy, nextx, cy + 50, n, bwidth, bheightS, gap);

    nextxi += 10; nextx += n * (bwidth + gap) + gap; n = 3;
    draw_chain(g, nextxi, cy, nextx, cy + 50, n, bwidth, bheightL, gap);
  }

  function draw_phase3(svg, sx, sy, width, height, bwidth, bheightL, bheightS) {
    var g = svg.append('g').attr('transform', 'translate(' + sx + ',' + sy + ')');
    var cx = width / 2, cy = 30;
    // g.append('rect').attr({ width: width, height: height, fill: 'lightgray', });
    g.append('polygon').attr({
      points: cx + ',5 ' + (cx + 250) + ',45 ' + (cx - 250) + ',45',
      style: "fill:blue; stroke:black; stroke-width:1"
    });

    var nextxi = cx - 200, nextx = cx - width / 2 + 20, n = 1, gap = 15;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap, true);

    for (var i = 0; i < 5; i++) {
      nextxi += 10; nextx += n * (bwidth + gap) + gap / 2; n = 0;
      draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap);
    }

    nextxi += 10; nextx += n * (bwidth + gap) + gap; n = 1;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightL, gap);

    for (var i = 0; i < 3; i++) {
      nextxi += 10; nextx += n * (bwidth + gap) + gap / 2; n = 0;
      draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap);
    }

    nextxi += 10; nextx += n * (bwidth + gap); n = 1;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap, true);

    nextxi += 10; nextx += n * (bwidth + gap); n = 1;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap);

    for (var i = 0; i < 4; i++) {
      nextxi += 10; nextx += n * (bwidth + gap) + gap / 2; n = 0;
      draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap);
    }

    nextxi += 10; nextx += n * (bwidth + gap); n = 1;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap, true);

    nextxi += 10; nextx += n * (bwidth + gap); n = 1;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap);

    for (var i = 0; i < 10; i++) {
      nextxi += 10; nextx += n * (bwidth + gap) + gap / 2; n = 0;
      draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap);
    }

    nextxi += 10; nextx += n * (bwidth + gap); n = 1;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap);

    for (var i = 0; i < 4; i++) {
      nextxi += 10; nextx += n * (bwidth + gap); n = 0;
      draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap);
    }

    nextxi += 10; nextx += n * (bwidth + gap) + gap; n = 2;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightL, gap);

    for (var i = 0; i < 4; i++) {
      nextxi += 10; nextx += n * (bwidth + gap); n = 0;
      draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap);
    }

    nextxi += 10; nextx += n * (bwidth + gap); n = 1;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightS, gap, true);

    nextxi += 10; nextx += n * (bwidth + gap); n = 1;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightL, gap, true);

    nextxi += 10; nextx += n * (bwidth + gap); n = 0;
    draw_chain(g, nextxi, cy + 10, nextx, cy + 50, n, bwidth, bheightL, gap, true);
  }

  function draw_chain(svg, cx, cy, sx, sy, n, width, height, gap, sorted) {
    if (n == 0) {
      svg.append('circle').attr({ cx: cx, cy: cy, r: 3, fill: 'white', stroke:'black' });
      return;
    }

    svg.append("path").attr({
      d: d3.svg.line().x(x).y(y).interpolate("basis")([
          { "x": cx,  "y": cy },
          { "x": cx,  "y": cy + 20 },
          { "x": sx,  "y": sy - 30 },
          { "x": sx,  "y": sy + width / 2 },
        ]),
      stroke: 'black',
      'stroke-width': 1,
      fill: 'none',
    });

    svg.append('circle').attr({ cx: cx, cy: cy, r: 3, fill: 'white', stroke:'black' });

    var g = svg.append('g');
    sx -= width / 2;
    for (var i = 1; i < n; i++, sx += width + gap) {
      draw_bucket(g, sx, sy, width, height, 100, sorted);
      draw_leaf_connectors(g, sx + width / 2, sy + height, sx + width * 3 / 2 + gap, sy);
    }
    draw_bucket(g, sx, sy, width, height, 50 + 50 * Math.random(), sorted);
  }

  function draw_bucket(svg, sx, sy, width, height, pfull, sorted) {
    svg.append('rect').attr({
      transform: 'translate(' + sx + ',' + sy + ')',
      width: width,
      height: height,
      fill: 'white',
      stroke: 'black',
      'stroke-width': 1,
    });

    svg.append('rect').attr({
      transform: 'translate(' + sx + ',' + sy + ')',
      width: width,
      height: height * pfull / 100,
      fill: sorted ? (height > 30 ? 'cyan' : 'blue') : 'gray',
    });
  }

  function draw_leaf_connectors(svg, x1, y1, x2, y2) {
    var dx = 15, dy = 15;
    svg.append("path").attr({
      d: d3.svg.line().x(x).y(y).interpolate("basis")([
          { "x": x1,  "y": y1 },
          { "x": x1,  "y": y1 + dy },
          { "x": x1 + dx,  "y": y1 + dy },
          { "x": x2 - dx,  "y": y2 - dy },
          { "x": x2,  "y": y2 - dy },
          { "x": x2,  "y": y2 },
        ]),
      stroke: 'black',
      'stroke-width': 1,
      fill: 'none',
    });
  }

  function x(d) { return d.x; }
  function y(d) { return d.y; }
})();
