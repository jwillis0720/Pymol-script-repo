#!/usr/bin/env python
def main():
	import urllib2
	import json
	api_key = 'VtxgIC2UnhfUmXe_pBksov7-lguAQMZD'
	url = 'http://www.energyhive.com/mobile_proxy/getCurrentValuesSummary?token='+api_key
	response = urllib2.urlopen(url)
	content = response.read()
	for x in json.loads(content):
	    if x["cid"] == "PWER":
	        print (x["data"][0].values()[0])


if __name__ == '__main__':
	main()

