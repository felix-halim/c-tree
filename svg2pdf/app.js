var express = require('express');
var app = express();
var fs = require('fs');
var path = require('path');
var exec = require('child_process').exec;

app.use(express.bodyParser());

app.get('/hello.txt', function(req, res){
  var body = 'Hello World';
  res.setHeader('Content-Type', 'text/plain');
  res.setHeader('Content-Length', body.length);
  res.end(body);
});

app.post('/p', function(request, response){
  console.log(request.body.user.name);
  console.log(request.body.user.email);
});

app.post('/svg2pdf', function(req, res) {
  res.setHeader('Access-Control-Allow-Headers', 'X-Requested-With');
  res.setHeader('Access-Control-Allow-Methods', 'GET, POST, OPTIONS');
  res.setHeader('Access-Control-Allow-Origin', '*');

  var prefix = '/tmp/' + (Math.random() * 1e20);
  console.log(prefix);

  fs.writeFile(prefix + '.svg', req.body.p, function (err) {
    if (err) {
      console.log(err);
      res.download(__dirname + '/test.pdf');
    } else {
      console.log("Converting to pdf");
      exec("~/bin/svg2pdf " + prefix + '.svg ' + prefix + '.pdf', function (error, stdout, stderr) {
        console.log("done");
				if (fs.existsSync(prefix + '.pdf')) {
	        res.download(prefix + '.pdf');
				} else {
	        res.download(prefix + '.svg');
				}
      });
    }
  }); 
});

app.listen(8000);
console.log('Listening on port 8000');
