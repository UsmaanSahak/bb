$(document).keyup(function() {
 if (document.activeElement.id == "searchLocus") {
  $("#CompleteSym").hide();
  $("#CompleteLocus").show();
  var query = document.getElementById("searchLocus").value;
  getSuggestions(query,"loSuggestions");
 } else if (document.activeElement.id == "searchSym") {
  var query = document.getElementById("searchSym").value;
  getSuggestions(query,"symSuggestions");
  $("#CompleteLocus").hide();
  $("#CompleteSym").show();
 } else {
  $("#CompleteLocus").hide();
  $("#CompleteSym").hide();
 } 
});


/*Make a jquery function on clicking a table entry (name it a class), in which it passes function
the value of the element that was clicked, and the function is what ends up calling getparalog().*/


$(document).click(function() {
 if (document.activeElement.id != "searchLocus" && document.activeElement.id != "searchSym") {
  $("#CompleteLocus").hide();
  $("#CompleteSym").hide();
 }
});


function getSuggestions(e,name) {
 var xhttp = new XMLHttpRequest();
 xhttp.onreadystatechange = function() {
  if (this.readyState == 4 && this.status == 200) {
   /*document.getElementById("loSuggestions").innerHTML = name;*/
   document.getElementById(name).innerHTML = this.responseText;
  }
 };
 var Url = "../cgi-bin/a.pl?entry=" + e + "&suggType=" + name;
 xhttp.open("GET",Url,true);
 xhttp.send();
}
































var gids, nG1=0, nG2, strainData, taxoData, dbLink, pubData, pfamData, mainRoleData, roleData, tigrData, acc_genome={}, acc_pls={}, groupName, orfData, geoData, sym_cid, treeData, gidC, gidP, uvIDs, conFuse=[], cid_order, cid_class, group_id;

var stdOid, orfPara, pathFun;

var r=25, startMax, endMax, stdStrain, zeroPoint;
var hUnit=16, hUnitL=16, hUnitR=13, hUnitP=23.6, tree_wCore=140, tree_wPls=105, wsvgPls=1024-tree_wPls,
 bgStartL=27, yEdge=28, base_tick=yEdge-12, hTick=4, unitW=6.5, lowerYbound=6, gapCont=60;
var hORF=5, lTri=8, hRNA=1.5;
var strand, guideline, guidelineP, ratioP, orfTranslate=[];

var spH=16, wTreeP=170, wThumbnail=100, color_paraName='#d9fdd9';

var transData, groupT, rangeT, wBar=80, hBlock=4, hUnitT=20, topEdge=26, guidelineT, transTitle;

$(document).ready(function(){
    $("#main_tabs").tabs();
    
    $.ajax({
        async:false, dataType:"json", url:"js-css/orf.json",
        success: function(data) { orfData = data }
    });

    $.ajax({
        async:false, dataType:"json", url:"js-css/data.json",
        success: function(data) {
            gids=data.gids;
            geoData=data.geo;
            taxoData = data.taxo;
            strainData=data.strain;
            dbLink=data.dbLink;
            pubData=data.pub;
            sym_cid = data.sym;

            for (var i=0; i<gids.length; i++){ if (gids[i]<1000){nG1++} }
            nG2 = gids.length - nG1;
            strainTbl();

            groupName = data.groupName;
            group_id = Object.keys(groupName);

            // get acc_genome
            $.each(orfData, function(acc,obj){
                var geno = obj.ID[0],
                    grp = obj.ID[1];
                if (grp[1]) {
                    conFuse.push(acc);
                    for (var i=0; i<grp.length; i++){ get_acc_g(geno,grp[i],acc) }
                } else {
                    get_acc_g(geno,grp,acc)
                }
            });
            
            gidC = data.gidC;
            gidP = data.gidP;
            transTitle=data.transTitle;
            transData = data.trans;
            groupT = data.group_locus;
            rangeT = data.range
        }
    });

    $.ajax({
        async: false, dataType: "json", url: "js-css/tree.json",
        success: function(data) {
            treeData = data;
 			draw_tree('#ltree')
		}
    });

    draw_tree('#gtree');
    repliconTbl();

    $.ajax({
        async:false, dataType:"json", url:"js-css/cdhit.json",
        success: function(data) {
            cid_order = data.cid_order;
            uvIDs = data.uv
        }
    });

// set orf colors by cdhit:
    cid_class = {};
    for (var i=0; i<cid_order.length; i++) { var j = (i+17)%17; cid_class[cid_order[i]] = 'cid' + (j+1) }

// pathFun
    var ypath=[];
    for (var i=0; i<gidC.length; i++){ ypath.push(yEdge + 9 + hUnit*i + (i? (i==ypath.length-1? 6 : 0) : (-6))) }
    var y1 = ypath.map(function(x){return x}).reverse();    
    $.merge(ypath,y1);
    pathFun = d3.svg.line()
        .x(function(d){return d})
        .y(function(d,i){return ypath[i]})
        .interpolate("linear");


    loadContigMenu();
    $.ajax({
        async: false, dataType: "json", url: "js-css/treeC.json",
        success: function(data) {draw_ctree(data)}
    });
    drawORF();

    // get acc_pls
    var tmpobj = {};
    $.each(orfData, function(acc,obj){
        var geno = obj.ID[0];
        if (geno>1000){ return;}
        var grp = obj.ID[1];
        if (grp==1 || grp==2 || grp==3){ return }
        var ac_gp = {};

        if (tmpobj[geno]){ ac_gp = tmpobj[geno]};
    
        if (grp[1]){ grp = grp[0] }
        ac_gp[acc] = grp;
        tmpobj[geno] = ac_gp
    });
    for (var i=0; i<gidP.length; i++){
        var gid = gidP[i];
        var ac_gp = tmpobj[gid];
        var aa = Object.keys(ac_gp);
        acc_pls[gid] = aa.sort(function(a,b){ return ac_gp[a] - ac_gp[b] })
    }

    $.ajax({
        async: false, dataType: "json", url: "js-css/treeP.json",
        success: function(data) {draw_ptree(data)}
    });
    drawPlasmid();

    loadRepliconMenu();
    addTransTitle();
    addTransContainer();
    drawSigLeg();
    draw_transbar(2);

    $("#sub_tabs").tabs();
    $('#treeP').css("width", wTreeP+'px');
    $('#plasmidZone').css("width", (wsvgPls-4)+'px');
    $('#paraNameC').css("margin-top", (spH*2.5+3)+'px');
    $('#treesC').css("margin-top", (spH*2.5)+'px');
	
	$('#replicon').change(function(){ drawORF() });
	drawLeg();

	$('#repliconMenu').change(function(){ draw_transbar($('#repliconMenu').val()) });
    
	$('#searchLocus').bind('keypress', function(e) { var code = e.keyCode || e.which; if(code == 13) { getParalog() } });
	$('#searchSym')  .bind('keypress', function(e) { var code = e.keyCode || e.which; if(code == 13) { getParalog(1)} });

	$("#seqP").mousemove(function(e) {
		var xMouse = e.pageX - this.offsetLeft - $('#nameP').width() + $("#seqP").scrollLeft()-190;
		var yMouse = e.pageY - this.offsetTop - 151;

	   	var ymax = $("#seqParea").height()-15;
    	$("#gdHp, #gdHp2").css("top", (yMouse>0? (yMouse<ymax? yMouse : ymax) : 0) +"px");
    	$("#gdHp").css("width", $("#seqP svg").width()+"px");
    	$("#gdHp2").css("width", $("#paraName").width()+"px");

	    if (xMouse >= $("#seqP svg").width()-5 || xMouse < 0) {return}
    	$("#gdVp").css("left", xMouse+"px");
        
	    var px = parseInt(xMouse/unitW+1);
	    guidelineP.attr("x1",px/ratioP).attr("x2",px/ratioP)
	});

    
    $("#seqO").mousemove(function(e) {
		var xMouse = e.pageX - this.offsetLeft - 15 + $("#seqO").scrollLeft();
		var yMouse = e.pageY - this.offsetTop - 102;

	   	var ymax = $("#seqOrth").height()+16;
    	$("#gdH").css("top", (yMouse>30? (yMouse<ymax? yMouse : ymax) : 30) +"px").css("width", $("#seqOrth svg").width()+"px");
    	$("#gdH2").css("top", (yMouse>30? (yMouse<ymax? yMouse-27 : ymax-27) : 3) +"px").css("width", $("#trees").width()+"px");

	    if (xMouse >= $("#seqOrth svg").width()-5 || xMouse < 0) {return}
    	$("#gdV").css("left", xMouse+"px");

    	var p = parseInt(xMouse/unitW+1);
    	$("#gMark").css("left", xMouse+"px").html(numberWithCommas(p));

	    var px = strand*parseInt(xMouse/unitW+1) + zeroPoint;
	    guideline.attr("x1",px/r).attr("x2",px/r)
	});

    $("#svgTrans").mousemove(function(e) {
		var yMouse = e.pageY - this.offsetTop - 147;
        if (yMouse<0 || yMouse > $("#guidelineT_area").height()-70){ yMouse=-50 }
	    guidelineT.attr("y1",yMouse).attr("y2",yMouse)
	});

    $('#close').on("click", function(){
    	$('#fasta').hide();
		d3.selectAll('#strain_tbl a').classed('act', false);
    	return false
    });

    $('#seqWin span').on("click", function(){
    	$('#seqWin').hide();
    });
    
    $("#main_tabs").tabs({
        activate: function(event, ui) {
            var activeTab = $('#main_tabs').tabs('option', 'active');
            if (activeTab==2 || activeTab==3){ $('#legClick').show() }
            else { $('#legClick').hide() }
        }
    });
    
    $(".down").menu();
    initMap();

    $.ajax({
        async:false, dataType:"json", url:"js-css/fam.json",
        success: function(data) {
            pfamData=data.pfam;
            mainRoleData=data.mainRole;
            roleData=data.role;
            tigrData=data.tigr
        }
    });

});


function get_acc_g(g,p,a){
    var obj = {};
    if (acc_genome[g]){ obj = acc_genome[g] }
    accs = [];
    if(obj[p]){ accs = obj[p] }
    accs.push(a);
    obj[p] = accs;
    acc_genome[g] = obj
}

function loadContigMenu(){
    for (var i=1; i<=3; i++){
        $('#replicon').append('<option value="'+i+'"' + (i==2? '  selected="selected"' : '') + '>'+ groupName[i][0] + '</option>')
    }
}

var matchAc={}, selectAc, repStart={}, svgP, rectBg, xscaleP, xscale0, lMax=0, lM=66242, ws=wsvgPls-4, point1={}, point2={}, base0=20, uvPls;
function drawPlasmid(){
    var gap=680, repH=7, orfH=6, gapZoomed=240, lAxis=0;
    for (var i=0; i<gidP.length; i++){
        var gid = gidP[i];
        var accs = acc_pls[gid];

        var acInGp = {};
        for (var j=0; j<accs.length; j++){
            var ac = accs[j];
            var gpid = orfData[ac].ID[1];
            gpid = gpid[1]? gpid : [gpid];
            for (var k=0; k<gpid.length; k++){
                var gp = gpid[k];
                var aa = [];
                if (acInGp[gp]){ aa=acInGp[gp] }
                aa.push(ac);
                acInGp[gp] = aa
            }
        }
        
        $.each(acInGp, function(gp,aa){
            if (aa.length==1){ return }
            var start=0;
            for (var j=0; j<aa.length; j++){
                var ac = aa[j];
                if (j){ repStart[ac] = start }
                start += orfData[ac].ID[2] + gapZoomed
            }
        });

        var l_gp = Math.max.apply(null, accs.map(function(d){return orfData[d].ID[2] + (repStart[d]? repStart[d] : 0)}));
        if (l_gp>lAxis){lAxis=l_gp}
//        if (l_gp>lM && gp!=33){lM=l_gp}

        var l=0;
        for (var j=0; j<accs.length; j++){
            var l0 = orfData[accs[j]].ID[2];
            l += l0
        }
        l += (accs.length-1)*gap;
        if (l>lMax){lMax=l}
    }
    
    svgP = d3.select("#pls").append("svg")
            .attr("width",ws).attr("height",hUnitP*gidP.length + base0 + 4);
    rectBg = svgP.append("rect")
            .attr("y",2)
            .attr("height",hUnitP*gidP.length+base0-2)
            .attr("id", "plsBg").attr("class", "hidden");
    xscaleP = d3.scale.linear().domain([0, lM]).range([0, ws]);
    xscale0= d3.scale.linear().domain([0, lMax]).range([0, wsvgPls]);

    tick_mark = svgP.append("g");
    drawAxis(tick_mark, lM/ws, lAxis, 1000, 5, 1000);
    svgP.append("rect").attr("id","cover").attr("class","whiteFill").attr("width",ws).attr("height",base0);

    var tick_mark = svgP.append("g").attr("id", "axisRep");
    drawAxis(tick_mark, lMax/wsvgPls, lMax, 10000, 5, 1000);

    uvPls = svgP.append("g").attr("class", "univShd");
    
    var base = hUnitP/2+base0;

    for (var i=0; i<gidP.length; i++){
        var gid = gidP[i];
        var accs = acc_pls[gid];
        
        var start0=0;
        $.each(accs, function(j,ac){
            var acData = orfData[ac];
                            
            var l = acData.ID[2],
                gpid = acData.ID[1];
 
            var tiptxt = gpid[1]? $.map(gpid, function(c){return groupName[c][0]}).join('|') : groupName[gpid][0];
            var plasmidRep = svgP.append("g")
                    .attr("class", 'plasmidRep')
                    .attr("id", ac)
                    .attr("transform", "translate(" + xscale0(start0) + "," +(base-repH/2) + ")")
                    .attr("start0",start0)
                    .attr("y0", base-repH/2)
				    .on("mouseover", function() {
                        $("#plstip").css("left", (d3.event.pageX-20)+"px").css("top", (d3.event.pageY-75)+"px").show().text(tiptxt);
                    })
                    .on("mouseout", function() {
                        $("#plstip").hide();
                    });
            
            if (!gpid[1]){ gpid = [gpid]}
            var fuseN = gpid.length,
                y0 = 0;
            for (var n=0; n<fuseN; n++){
                plasmidRep.append("rect")
                    .attr("class", "fillLess orth"+(gpid[n]==33? 32 : gpid[n]-3))
                    .attr("width", xscale0(l)).attr("height", repH/fuseN)
                    .attr("transform", "translate(0," + y0 + ")")
                    .on("mouseover", function(){
                        var gp = $(this).attr("class").split('orth')[1].split(' ')[0];
                        d3.selectAll('.orth'+gp).classed("fillLess",false)
                    })
                    .on("mouseout", function(){d3.selectAll($(".plasmidRep rect")).classed("fillLess", true)})
                    .on("click", function(){
                        rectBg.classed("hidden",false);
                        selectAc=$.map($("."+$(this).attr("class")).closest("g"), function(a){return a.id});
                        msClick(0,$(this).attr("class").split('orth')[1].split(' ')[0])
                    });
                y0 += repH/fuseN
            }

            var orfStart = repStart[ac]? repStart[ac] : 0;
            var orfs = svgP.append("g")
                .attr("class","orfs hidden")
                .attr("id", 'o'+ac)
                .attr("transform", "translate(" + xscaleP(orfStart) + "," +(base-orfH/2) + ")");
            
            $.each(acData, function(loc,obj){
                if (loc=='ID'){ return }
                if (!obj.cid || obj.cid<0) { return }
                var p1 = (obj.L>0? obj.end-obj.L : obj.end),
                    s1 = xscaleP(p1),
                    l1 = xscaleP(Math.abs(obj.L));
                var ids;
                if (obj.orth){
                    var ids = obj.cid+'_'+obj.orth;
                    var pp1=[], pp2=[], cons=[];
                    if (point1[ids]){ pp1 = point1[ids]; pp2 = point2[ids]; cons=matchAc[ids] }
                    cons.push(ac);

                    var s2 = s1 + l1;       
                    pp1.push([s1+xscaleP(orfStart), base-repH/2-orfH], [s1+xscaleP(orfStart), base+repH/2+orfH]);
                    pp2.unshift([s2+xscaleP(orfStart), base+repH/2+orfH], [s2+xscaleP(orfStart), base-repH/2-orfH]);
                    
                    point1[ids] = pp1;
                    point2[ids] = pp2;
                    matchAc[ids] = cons
                }
                
                var isUv;
                for (var n=0; n<fuseN; n++){
                    var gp = gpid[n];
                    if ($.inArray(ids,uvIDs[gp])!=-1){ isUv=1; break }
                }
                
                var lo = orfs.append("rect")
                    .attr("class", isUv? 'cidUniv' : cid_class[obj.cid])
                    .attr("transform", "translate(" + s1 + "," + ((repH+orfH)/2*(obj.L>0? 1 : -1)) + ")")
                    .attr("width", l1).attr("height", orfH)
                    .on("mouseover", function(){
                        var tipText = loc;
                        if (obj.sym) { tipText += '&nbsp;&nbsp;[ <em>' + obj.sym + '</em> ]' }
                        $("#plstip").css("left", (d3.event.pageX-32)+"px").css("top", (d3.event.pageY-75)+"px").html(tipText)
                            .show()//.fadeIn().delay(400).fadeOut()
                    })
                    .on("mouseout", function(){$("#plstip").hide()})
                    .on("dblclick", function(){getParalog(2,obj.cid,loc)});

                if (ids){ lo.on("click", function(){msClick(ids)}) }
                
                if (!obj.orth){ lo.style("fill", "white") }  
            });

            start0 += l + gap
        });
        base += hUnitP
    }
    var rib = svgP.append("polygon").attr("id","rib").attr("class", "shade hidden")        
}

function msClick(ids,rep){
    d3.selectAll('.univShd polygon').classed("hidden", true);
    if (ids){selectAc = matchAc[ids]}
    
    var maxW = Math.max.apply(null, selectAc.map(function(id){return orfData[id].ID[2] + (repStart[id]? repStart[id] : 0)}));
    maxW = xscaleP(maxW) + 4;
    rectBg.attr("width", maxW-2);
            
    d3.selectAll(".plasmidRep")
        .transition()
            .duration(1500)
            .attr("transform", function(){
                var myy = $(this).attr("y0");
                if ($.inArray(this.id,selectAc)!=-1){
                    return "translate(" + (repStart[this.id]? xscaleP(repStart[this.id]) : 0) + "," + myy + "), scale(" + lMax/lM + ",1)"
                }
                var myx = xscale0($(this).attr("start0"));
                return "translate(" + (myx+maxW) + "," + myy + ")"
            });
    svgP.attr("width",ws+maxW);
    d3.selectAll("#axisRep, #cover").transition().duration(1500).attr("transform", "translate(" + maxW + ",0)");

    if (rep){
        rep = rep*1+3;
        if (uvIDs[rep]){
            var uvco=uvIDs[rep];
            if (!d3.selectAll('.uvpls'+rep)[0][0]){
                for (var i=0; i<uvco.length; i++){
                    var id = uvco[i],
                        points1 = $.map(point1[id], function(x,i){ return x.join(',') }).join(' '),
                        points2 = $.map(point2[id], function(x,i){ return x.join(',') }).join(' ');
                    uvPls.append("polygon").attr("class","uvpls" + rep)
                        .attr("points", points1 + ' ' + points2);
                }
            }
        }
    }
    
    d3.selectAll(".orfs")
        .classed("hidden", function(){ return $.inArray(this.id.substring(1),selectAc)==-1 })
        .transition()
            .delay(1000)
            .duration(1200)
            .style('opacity', function(){return $.inArray(this.id.substring(1),selectAc)==-1? 0:1})

    if (rep){
        d3.selectAll(".uvpls"+rep)
            .classed("hidden", false).style('opacity',0)
            .transition().delay(1000).duration(1200).style('opacity',1)
    }

    d3.select('#rib').classed("hidden", true).style('opacity', 0);
    
    if (ids){
        var points1 = $.map(point1[ids], function(x,i){ return x.join(',') }).join(' '),
            points2 = $.map(point2[ids], function(x,i){ return x.join(',') }).join(' ');
        d3.select("#rib").attr("points", points1 + ' ' + points2)
            .classed("hidden", false)
            .transition().delay(600).duration(1200).style('opacity', 1)
    }
}


function locateORF(cid){
    d3.selectAll('#syn svg .shade:not(#shadeC)').remove();
    
    for (var n=0; n<gidC.length; n++){
        var base = yEdge + 9 + hUnit*n,
            myorfs = $('#syn .' + n + '_' + cid);
        
        if (!myorfs[0]){ continue }
        
        $.each(myorfs, function(i,myorf){
            var points = myorf.attributes.points.value.split(' ');
            var p1 = points[0].split(',')[0],
                p2 = points[2].split(',')[0],
                px = Math.min(p1,p2);
            
            px += orfTranslate[n];
            
            var shd = d3.select('#syn svg').append("rect")
                .attr("class", 'shade dotStroke')
                .attr("x", px).attr("y", base-7)
                .attr("width", Math.abs(p2-p1)).attr("height", 14);
            
            if (!stdOid || stdOid && cid!=stdOid.split('_')[0]){
                shd.style("stroke", "gray")
            }
        })
    }
}

var cdhit, locus, symbol, idsP, locus_old, paraSeq, paraDnd, orthSeq;
function getParalog(searchSym, myId, myLocus) {
    $('#newsPara').css("left", (searchSym? 788 : 335) + 'px').hide();
    if (myId || (!myId && myLocus)){
        $('#main_tabs').tabs({ active:1 });

        if (myLocus) {
            if (locus && locus==myLocus) { $('#searchLocus').val(locus); $('#searchSym').val(symbol); return }
            locus = myLocus
        }
        if (cdhit && myId && cdhit==myId) {
            if (myLocus){ changeHL() }
            $('#searchLocus').val(locus); $('#searchSym').val(symbol);
            return
        }
        if (myId){ cdhit = myId }
    } else {
        var v = $(searchSym? '#searchSym' : '#searchLocus').val();
        v = v.replace(/ /g, '');
    
        if (v==''){ $('#newsPara').html('Please input ' + (searchSym? 'gene symbol.' : 'locus name.')).show(); return }
        else { $('#newsPara').hide() }

        v = v.toUpperCase();

        if (!searchSym){
    	   if (locus && locus.toUpperCase()==v) { $('#searchLocus').val(locus); $('#searchSym').val(symbol); return }
        } else {
            if (v.indexOf(',')>0){ v = v.split(',')[0] }
            var cid, syms=[];
    	   $.each(sym_cid, function(s,c){
    		  if (s.toUpperCase() == v){ cid = c; return false }
    		  else if (v==s.toUpperCase().substring(0, v.length)) { syms.push(s) }
    	   });
    	   if (!cid) {
    		  var inf;
    		  if (syms.length==0) { inf = '<b>' + $('#searchSym').val() + '</b> not found in database' }
    		  else { inf = 'Please select: &nbsp;' + syms.sort().join('&nbsp;&nbsp; ') }
    		  $('#newsPara').html(inf).show();
    		  return
    	   }
    	   if (cdhit && cdhit==cid) { $('#searchLocus').val(locus); $('#searchSym').val(symbol); return }
    	   cdhit = cid
        }
    }

    $.ajax({
		url: "cgi-bin/get-paralog.pl",
        dataType: "json",
        data: 'ip=' + (searchSym? cdhit : (myLocus? locus.toUpperCase() : v)),
		error: function(){ alert('AJAX error.') },

		success: function(data) {
			if (!searchSym && !data.cdhit){
				var inf = data.noFound? 'not found in database' : (data.seq_err? 'has sequence error' : (data.rna? 'is ' + data.rna + 'RNA' : 'is pseudogene'));
                $('#newsPara').html('<b>' + (myLocus? myLocus : $('#searchLocus').val()) + '</b> ' + inf).show();
                return
            }
            $('#dbl_click').show();
            if (!myLocus){ locus = data.locus? data.locus : data.ids[0] }
            $('#searchLocus').val(locus);

            if (cdhit==data.cdhit){ changeHL(); $('#searchSym').val(symbol); return }
            if (!searchSym){ cdhit = data.cdhit }

            symbol = data.symbol? data.symbol.join(', ') : '';
			$('#searchSym').val(symbol);
            
			locus_old = locus.toUpperCase();

			$('#paraName tr').remove();
			$('#treeP').empty();

			idsP = data.ids;
			orfPara = data.orf;
			paraName();

            $('#gdHp2').hide();
            $('#paraInf').show();
			paraPfam(data.pfam);
			if (data.tigr) { paraTigr(data.tigr) } else { $('#tigr').hide() }
			if (data.pmed) { paraPub(data.pmed) } else { $('#pubmed').hide() }
            
            var wSubTabP = 1012-wTreeP-$('#nameP').width();
            $('#sub_tabs').tabs({ active:0 })
                          .css({'width': wSubTabP + 'px',
                                'min-height': (idsP.length*spH+39) + 'px'});

            if (cdhit<0){ $('#sub2').hide(); return }            
            $('#sub2').show();
            
            if (data.nodes){
                $('#downPara ul li').eq(2).show();
                draw_treeP(data.nodes, data.w_tree, data.unit)
            } else {
                $('#downPara ul li').eq(2).hide()
            }

            $('#seqP').css("width", (wSubTabP-wThumbnail-18) + 'px');
            
            paraSeq = data.seq;
            paraDnd = data.dnd;
        }
    });
    
    function changeHL(){
        var lu = locus.toUpperCase();
        if (lu!==locus_old){
            $('#' + locus_old + ', #paraNameC ' + locus_old).css("background-color", "inherit");
            $('#' + lu + ', #paraNameC ' + lu).css("background-color", color_paraName);
            locus_old = lu
        }
    }
}

function paraPfam(pfam){
    $('#pfam span').html(pfam? pfam : 'unclassified');
    $('#pfam tr:not(:first-child)').remove();
    $('#pfam th').css("padding-bottom", pfam && pfamData[pfam]? '2px' : '5px')
    if (pfam && pfamData[pfam]) { $('#pfam').append('<tr><td>' + pfamData[pfam] + '</td></tr>') }
}

function paraTigr(tigr){
	var trs;
    var n=0;
    $.each(tigr, function(i,id){
            var obj = tigrData[id];
            if (n){ trs += '<tr><td colspan="2"><hr></td></tr>' }
			trs += '<tr><td>TIGR ID</td><td>' + id + '</td></tr>';
			trs += '<tr><td>Description</td><td>' + obj[0] + '</td></tr>';
			if (obj[1]){ trs += '<tr><td>Comment</td><td>' + obj[1] + '</td></tr>' }
			if (obj[2]){ trs += '<tr><td>DB Ref</td><td>' + obj[2].join(' | ') + '</td></tr>' }
			if (obj[3]){ trs += '<tr><td>Pub Ref</td><td>' + obj[3].join(' | ') + '</td></tr>' }
			if (obj[4]){ trs += '<tr><td>Role</td><td>' + $.map(obj[4], function(v,i){
                    var myRole = roleData[v];
                    return (myRole[0]? ('Main: ' + mainRoleData[myRole[0]] + '; ') : '') + 'Sub: ' + myRole[1]
                }).join('<br>') + '</td></tr>'
            }
			if (obj[5]){ trs += '<tr><td>GO</td><td>' + $.map(obj[5], function(v,i){return '<a href="' + dbLink.go + v + '" target="_blank">GO:'+v + '</a>'}).join(' | ') + '</td></tr>' }
            n++
		});
	
	$('#tigr tr:not(:first-child)').remove();
	$('#tigr').append(trs).show()
}


function paraPub(pub){
	var trs;
    $.each(pub, function(n,id){
        var obj = pubData[id];
        trs += '<tr><td>&bull;</td><td>' + (obj[0]? obj[0] + ' <em>et al</em>. ' : '[No authors listed]. ') + obj[3] + '. ' + obj[1] + ' <em>' + obj[2] + '</em>. <a href="' + dbLink.pubmed + id + '" target="_blank">PMID' + id + '</a></td></tr>';
    });
	$('#pubmed tr:not(:first-child)').remove();
	$('#pubmed').append(trs).show()
}


function paraName(){
	var trs;
	$.each(idsP, function(i,id){
        var arr = orfPara[id];
        var orth = arr.pop();

        trs += '<tr height="' + spH + 'px" id="' + id.toUpperCase() + '"' + (orth? (' class="orth' + orth + '"') : '') + (id.toUpperCase()==locus_old? (' style="background-color:' + color_paraName + '"') : '') + '><td>';
        if (orth){
            var rep = arr[1][0], rp;
            if (rep==1 || rep==2 || rep==3){ if ($.inArray(arr[0],gidC)!=-1){ rp=2 } }
            else { rp=3 }
            if (rp){ trs += '<a href="#" ondblclick="gotoalign(' + rep  + ',' + rp + ',' + orth + ')">'}
        }
        trs += id;
        if (orth){ trs += '</a>' }
        trs += '</td><td>';

        if (arr[0]==120){ arr[0]=100 }
        var strainObj = strainData[arr[0]];
        var ar = [strainObj.strain, taxoData[strainObj.tid]];
        ar.push($.map(arr[1], function(c){ return groupName[c][0] }).join('|'));
		trs += ar.join('</td><td>');
		trs += '</td></tr>';
	});
	$('#paraName').append(trs);
    $('#paraName a').click(function(e){ e.preventDefault() })
}

function gotoalign(rep, rp, orth){
    $('#main_tabs').tabs({ active:rp });
    var ids = cdhit + '_' + orth;
    if (rp==2){
        if (rep != $('#replicon').val()){ $('#replicon').val(rep).change() }
        alignOrth(ids)
    }
    else { msClick(ids) }
}

var cid_seq;
function drawAlign(orth) {
	var seq, orderList, unitH, seqContainer, seqArea,
        colorAA = ["G","S","T","Y","C","Q","N","K","R","H","D","E","-","X"];
    
	if (orth){
		if ($('#treesC').is(':empty')){ $('#trees').clone().appendTo('#treesC') }
		$('#seqO,#gdH2,#downOrth,#treesC').show();
        $('#seqOrth').empty();

		seq = orthSeq;
		orderList = gidC;
		unitH = hUnit;
		seqContainer = 'seqOrth';
        seqArea = 'seqO'
	} else {
		$('#seqP, #paraNameC,#gdHp2').show();
		$('#paraName').clone().appendTo('#paraNameC');

        seq = paraSeq;
		orderList = idsP;
		unitH = spH;
		seqContainer = 'seqParea';
        seqArea = 'seqP'
    }

    var seq_length = seq[stdStrain][0].length,
        svg_w = seq_length*unitW,
        svg_h = unitH*orderList.length*2 + 2.5*spH;
    

    var svg = d3.select('#'+seqContainer).append("svg").attr("width", svg_w).attr("height", svg_h);

    var codenBg = svg.append("g").attr("class", 'codenBg');
    var noSeq = svg.append("g").attr("class", 'grayFill');
    
    for (var i=0; i<seq_length/6; i++){
        codenBg.append("rect").attr("x",i*6*unitW).attr("y",0).attr("width",unitW*3).attr("height",svg_h)
    }

    var adj = orth? -6 : -5;
    var base = adj;
    drawSeq(0);

    base += 4;
    svg.append("rect").attr("x",0).attr("y",base).attr("width",svg_w).attr("height",spH*3-3).style("fill","white");
    
    var tick = svg.append("g").attr("class", "tick"),
    	mark = svg.append("g").attr("class", "mark");
    
    base +=spH/2;
    for (var i=1; i<seq_length/30; i++){
        var px = i*30*unitW-3,
            py = base;
        tick.append("line").attr("x1",px).attr("x2",px).attr("y1",py).attr("y2", py-spH/2);
        py += spH-3;
        mark.append("text").attr("x",px+3).attr("y",py).text(i*30);
        py += spH-3;
        mark.append("text").attr("x",px+3).attr("y",py).text(i*10);
        py += spH/2+3;
        tick.append("line").attr("x1",px).attr("x2",px).attr("y1",py).attr("y2", py-spH/2)
    }
    
    base += spH*2-1;
    drawSeq(1);
    
    $("#"+seqArea).scrollLeft(0);
    
    function drawSeq(idx){
        var txt = svg.append("g").attr("class", idx? 'seqAA' : 'seqNt');
        var ss = seq[stdStrain][idx];
        $.each(orderList, function(n,loc){
            base += unitH;
            
            if (orth && !seq[loc]){
                var block = noSeq.append("rect")
                            .attr("x",0).attr("y",base-unitH-adj)
                            .attr("width",svg_w).attr("height",unitH);
                return true
            }
            
            var cc = seq[loc][idx];
            if (loc != stdStrain){ cc = matchSeq(ss, cc, idx) }
            var text = txt.append("text").attr("x", idx? unitW : 0).attr("y", base);
            if (idx){
                var arr = cc.split('');
                for (var i=0; i<arr.length; i++){
                    var textAA = text.append("svg:tspan").html(arr[i]=='.'? '&nbsp;' : arr[i]);
                    if ($.inArray(arr[i], colorAA) != -1){ textAA.attr("class", 'Lewis_'+arr[i]) }
                }
            } else {
                text.html(cc);
                if (loc==stdStrain){ text.attr("class", "stdNt") }                
            }
        });
    }
}

function matchSeq(s,c,idx){
    var arrS = s.split(''),
        arrC = c.split('');
    return $.map(arrC, function(a,i){ return a==arrS[i]? (idx? '.' : '&nbsp;') : a }).join('')
}

function hideClone() { $('#paraNameC, #gdHp2').hide() }


function downSeq(i, orth){
	var content='';
	if (i==2){ content += paraDnd }
	else {
		if (orth){
			$.each(gids, function(n,id){
				if (!orthSeq[id]){ return }
				content += '>' + strainData[id].strain + '|' + orthSeq[id][2] + '<br>' + orthSeq[id][i] + '<br>'
			})
		} else {
			$.each(idsP, function(n,id){content += '>' + orfPara[id][0] + '|' + id + '<br>' + paraSeq[id][i] + '<br>'})
		}
	}
	$('#seqDown').html(content);
	$('#seqWin').show();
}

function loadRepliconMenu(){
    $.each(groupT, function(id,obj){
        var name = groupName[id][0];
        if (id==32){ name = groupName[26][0] + '|' + name }
        $('#repliconMenu').append('<option value="'+id+'"' + (id==2? '  selected="selected"' : '') + '>'+ name + '</option>')
    });
}

function draw_transbar(c) {
    var transOrf = groupT[c];

    $('#transLocus').empty();
    $.each(transOrf, function(n,l){
        $('#transLocus').append('<p><a href="#" ondblclick="getParalog(0,0,\'' + l  + '\')">' + l + '</a></p>')        
    });
        
    var svg_h = topEdge + hUnitT * (transOrf.length);
    var y1 = topEdge-6, y2 = svg_h-hUnitT+6;

    for (var i=0; i<transTitle.length; i++){
        var infl = transTitle[i],
            inflRange = infl.replace(/[0-9]/g, ''),
            rg = rangeT[inflRange];
        
        $('#'+transTitle[i]).empty();

        var svg = d3.select('#' + transTitle[i]).append("svg").attr("width", wBar).attr("height", svg_h);

        var scale = d3.scale.linear().domain([-rg,rg]).range([1, wBar-2]);
        
        svg.append("rect").attr("class", "bg_trans")
	       .attr("x", scale(-rg))
           .attr("y", y1)
           .attr("width",scale(rg*2))
           .attr("height", y2-y1);
    
        //draw axis
        var xAxis = d3.svg.axis().scale(scale).orient("top").ticks(5).tickSize(-svg_h, 0);
        svg.append("g").attr("class", "xaxis")
            .attr("transform", "translate(0,18)")
            .call(xAxis);

        svg.append("line").attr("class", 'zeroline')
            .attr("x1", scale(0)).attr("x2", scale(0))
            .attr("y1", y1).attr("y2", y2);

        $.each(transOrf, function(n, l) {
            if (!transData[l] || !transData[l][infl]){ return }
            var vv = transData[l][infl];
            var base = topEdge + hUnitT * n;
        
            var v = i<8? vv[0] : vv;
            svg.append("line")
                .attr("x1", scale(0)).attr("x2", scale(v))
                .attr("y1", base).attr("y2", base)
                .attr("class", vv[1] || i>7? "transbarSig" : "transbarNo")
        });
    }
    
    //draw guideline
    $("#guidelineT_area").empty();
    var w = $('#svgTrans').width()? $('#svgTrans').width() : 1020;
    var svg = d3.select("#guidelineT_area")
              .append("svg")
              .attr("width", w-85)
              .attr("height", ($('#svgTrans').height()? $('#svgTrans').height() : 557)-10);
    guidelineT = svg.append("line").attr("class","guideline")
    			.attr("x1",0).attr("x2",w)
    			.attr("y1","-50px").attr("y2","-50px")
}

function addTransTitle(){
    var transTitle0 = [['c-di-GMP',25987708], ['Growth Phases',27706236,3], ['Fitness in Carbohydrates',27279039,4], ['Life stages',25425211,3]];
    var tdAdd;
    for (var i=0; i<transTitle0.length; i++){
        tdAdd += i? '<td class="spacer"></td>' : 0;
        tdAdd += '<th' + (transTitle0[i][2]? ' colspan="' + transTitle0[i][2] + '"' : '') + '><a href="' + dbLink.pubmed + transTitle0[i][1] + '" target="_blank">' + transTitle0[i][0] + '</a></th>';
    }
    $('#transTitle0').append(tdAdd)
}

function addTransContainer(){
    var trAdd = '<tr><td id="transLocus"></td>';
    for (var i=0; i<transTitle.length; i++){
        trAdd += '<td id="' + transTitle[i] + '" class="transBar"></td>';
        if (!i || i==3 || i==7){ trAdd += '<td></td>' }
    }
    trAdd += '</tr>';
    $('#svgTrans').append(trAdd);
    
    var div = $('<div>').attr('id', 'guidelineT_area');
    $('#svgTrans').append(div);
}

var conOld;
function alignOrth(oid) {
    if (oid == stdOid && conOld==$('#replicon').val()) { scrollSyn(); return }
	stdOid = oid;

	$('#seqO, #gdH2, #downOrth, #treesC').hide();
    $('#btnShowSeqO').show()
 
  // length of fig, start point:
    startMax = 0;
    var remainMax = 0;
    $.each(gidC, function(n, gid) {
        var myid = '#'+gid+'_'+stdOid;
        if (!$(myid)[0]) { return }
        var myorf = $(myid)[0];
        var end = myorf.attributes['end'].value;

        if (end > startMax) { startMax = end*1; endMax = end - myorf.attributes['len'].value }
        var remain = $('.d'+gid)[0].attributes['len'].value - end;
        if (remain > remainMax) { remainMax = remain }
    });

    var lMax = startMax + remainMax;
    d3.select('#syn svg').attr("width", lMax/r+1);

    //prth
    var uvco=uvIDs[$('#replicon').val()], isUv, pathT={};
    if ($.inArray(stdOid, uvco)!=-1){ isUv=1 }

    stdStrain = '';
    var firstGid = '', shadeX = [], shadeY = [], ii=0;
    $.each(gidC, function(n, gid) {
    	var myid = '#'+gid+'_'+stdOid;
		if (!$(myid)[0]) { return }
		var myorf = $(myid)[0];
		var end =  myorf.attributes['end'].value;
		var len =  myorf.attributes['len'].value;

		var base = yEdge + 9 + hUnit * n;
        
        orfTranslate[n] = Math.round((startMax-end)/r*10)/10;
        d3.select('.d'+gid)
            .attr("transform", "translate(" + orfTranslate[n] + ")");
        if (!firstGid) {
            d3.select("#tick_mark")
                .attr("transform", "translate(" + orfTranslate[n] + ")");
            strand = len>0? 1 : -1;
            firstGid = gid
        }

        // path
        if (isUv){
            for (var i=0; i<uvco.length; i++){
                var id = uvco[i],
                    arr = xpath[id];
                var arrNew=[];
                if (pathT[id]){ arrNew = pathT[id] };
                arrNew[ii] = arr[ii].map(function(x){ return Math.round((x+orfTranslate[n])*10)/10 });
                pathT[id] = arrNew
            }
        }
        ii++;
        
	    //  calculate shade points
	    if (!shadeX.length) { shadeX.push(startMax, (startMax-len)); shadeY.push(base-7, base-7) }
	    shadeX.push(startMax-len); shadeY.push(base);

	    // standard gid, zeroPoint
	    if (!stdStrain) { stdStrain = gid; zeroPoint = startMax - len }
    });

    //path
    d3.selectAll('.univShd path').classed("hidden", isUv? false : true);
    if (isUv){ uvShade(pathT) }
    
    var lastX = shadeX[shadeX.length-1], lastY = shadeY[shadeY.length-1];
    shadeX.push(lastX, shadeX[0]); shadeY.push(lastY+7, lastY+7);

    //  draw shade
    var points = $.map(shadeX, function(vx,i){ return vx/r + ',' + shadeY[i]}).join(' ');
    d3.selectAll('#syn svg .shade').remove();
    d3.select('#syn svg').append("polygon").attr("id", "shadeC").attr("class","shade").attr("points", points);

    scrollSyn();
    guideline.attr("x1","-50px").attr("x2","-50px");
    conOld = $('#replicon').val();
}


function scrollSyn() { $("#showSynteny").scrollLeft((startMax+endMax)/2/r-389/2) }

function getAlign() {
    $('#btnShowSeqO').hide();
	
    $.ajax({
		url: "cgi-bin/get-ortholog.pl",
        dataType: "json",
		data: "id=" + stdOid,
		error: function(){ alert("AJAX error.") },

		success: function(data) {
			orthSeq = data;
			$.each(gidC, function(n,gid) {
				if (orthSeq[gid]){ stdStrain=gid; return false }
			});
	    	drawAlign(1)
		}
    });
}

function thumbnail(){
    if (cid_seq && cid_seq==cdhit) { return }
    cid_seq = cdhit;

    $('#seqP').hide();
    $("#thumbnail, #seqParea").empty();
    
    stdStrain = idsP[0];
    var pLength = paraSeq[stdStrain][0].length;
    ratioP = pLength/wThumbnail;
    var svg = d3.select("#thumbnail") .append("svg") .attr("width", wThumbnail+1) .attr("height", lowerYbound+(idsP.length+2)*spH);

    var line = svg.append("g").attr("class", 'tnStroke');
    var seg = svg.append("g").attr("class", 'orfPara');
    
    $.each(idsP, function(i,id){
        var py = (spH-hORF)/2+spH*i+1;
        line.append("rect")
            .attr("x", 0)
            .attr("y", py)
            .attr("width", wThumbnail)
            .attr("height", hORF);
        
        var seq = paraSeq[id][0];
        var ss = [], good, blank;
        for (var n=0; n<seq.length; n++){
            if (seq.substring(n,n+1)!='-' && !good){ ss.push(n); good=1; blank=0 }
            else if (seq.substring(n,n+1)=='-' && !blank) {
                good=0; blank=1;
                if (n) { ss.push(n-1) }
            }
        }
        for (var j=0; j<ss.length; j+=2){
            var s1 = ss[j];
            var s2 = ss[j+1]? ss[j+1] : seq.length;
            seg.append("rect")
                .attr("x", s1/ratioP)
                .attr("y", py)
                .attr("width", (s2-s1)/ratioP)
                .attr("height", hORF)
        }
    });
    
    var base = (idsP.length+1)*spH;
    line.append("line").attr("x1", 0).attr("x2", wThumbnail).attr("y1", base).attr("y2", base);
    line.append("line").attr("x1", 0).attr("x2", 0).attr("y1", base-2).attr("y2", base+2);
    line.append("line").attr("x1", wThumbnail).attr("x2", wThumbnail).attr("y1", base-2).attr("y2", base+2);

    svg.append("text")
        .attr("x",wThumbnail/2)
        .attr("y", base+12)
        .text(numberWithCommas(pLength)+' bp');

//draw guideline
    guidelineP = svg.append("line").attr("class","guideline")
    			.attr("y1",0).attr("y2",base)
    			.attr("x1","-50px").attr("x2","-50px")
}

var xpath;
function drawORF() {
    var grp = $('#replicon').val(),
        uvco=uvIDs[grp];
    xpath = {};

	$('#btnShowSeqO, #seqO, #gdH2, #downOrth').hide();
    $('#treesC').empty();
    
  // length of fig:
    var lMax = 0;
    for (var n=0; n<gidC.length; n++){
        var gid = gidC[n],
            accs = acc_genome[gid],
            acc = accs[grp];
        var l = 0;
        $.each(acc, function(i,ac){
            var conPara = orfData[ac].ID;
            l += conPara[2];
            if (conPara[3]){ l += conPara[4]? conPara[4] : gapCont }
        })
        if (l>lMax){ lMax = l }
    }
    
    var svgH = hUnit*gidC.length+yEdge+5;
    
    $("#syn").empty();
    var svg = d3.select("#syn").append("svg").attr("width", lMax/r+1).attr("height", svgH),
        univShd = svg.append("g").attr("class", 'univShd'),
        xscale = d3.scale.linear().domain([0, 100]).range([0, 100/r]);

    for (var i=0; i<uvco.length; i++){
        univShd.append('path').attr("id", 's'+uvco[i])
    }
    
    
    $.each(gidC, function(n, gid) {
        orfTranslate[n] = 0;
        
		var base = yEdge + 9 + hUnit * n;
        
        var accs = acc_genome[gid],
            acc = accs[grp];
        if (acc.length>1){
            acc = acc.sort(function(a, b) { return (orfData[a].ID[3]? orfData[a].ID[3] : 0) - (orfData[b].ID[3]? orfData[b].ID[3] : 0) })
        }

		var myStrain = svg.append("g").attr("class", 'd'+gid);
        var line = myStrain.append("g").attr("class", "baseline"),
			tag = myStrain.append("svg:a").attr("xlink:href", "#");

        var start0 = 0;
        $.each(acc, function(i, ac){
            var oo = orfData[ac];
                        
            //draw base line
            var pp1 = start0;
            if (oo.ID[3]){ pp1 += oo.ID[4]? oo.ID[4] : gapCont }
            var pp2 = pp1 + oo.ID[2];
            
            var p1 = xscale(pp1), p2 = xscale(pp2);

            line.append("line").attr("x1", p1).attr("x2", p2).attr("y1", base).attr("y2", base);
	    	line.append("line").attr("x1", p1).attr("x2", p1).attr("y1", base-hORF*1.2).attr("y2", base+hORF*1.2);
	    	line.append("line").attr("x1", p2).attr("x2", p2).attr("y1", base-hORF*1.2).attr("y2", base+hORF*1.2);
            
 		//draw orf
           $.each(oo, function(id,obj){
                if (id=='ID'){ return }
                var end = Math.round(xscale(obj.end + pp1)*10)/10,
                    start = Math.round(xscale(obj.end-obj.L + pp1)*10)/10;
                var xv = [start, end-lTri*(obj.L>0? 1 : -1), end, end-lTri*(obj.L>0? 1 : -1), start],
                    hB = obj.rna? hRNA : hORF,
                    yv = [base-hB/2, base-hB/2, base, base+hB/2, base+hB/2],
                    points = $.map(xv, function(x,i){ return x + ',' + yv[i]}).join(' ');

                var tipText = id;
                if (obj.sym) { tipText += '&nbsp;&nbsp;[ <em>' + obj.sym + '</em> ]' }
               
               var co;
               if (obj.cid>0 && obj.orth){ co = obj.cid + '_' + obj.orth }
               
               var isUv = co && $.inArray(co, uvco)!=-1;
                var myclass =  obj.cid>0 ? (isUv? 'cidUniv' : cid_class[obj.cid]) : 'cidna';
                if (obj.cid && !obj.orth){ myclass += ' ' + n+'_'+obj.cid }

                var thisOrf = (obj.cid? tag : myStrain).append("polygon")
				    .attr("points", points)
				    .attr("class", myclass)
				    .on("mouseover", function() {
				        var xMouse = d3.event.pageX + $("#showSynteny").scrollLeft() - tree_wCore - 40;
				        $("#syntip").css("left", xMouse+"px").css("top", (base-23) + "px").show().html(tipText)
                    })
                    .on("mouseout", function() { $("#syntip").hide() });
               
                if (!obj.orth) {
                    thisOrf.style("fill", "white");
                    if (obj.pseu) { thisOrf.style("stroke-dasharray", "4,1") }
			     }

                if (!obj.cid){ return }
                
                thisOrf.on("dblclick", function(){ getParalog(2,obj.cid,id) });
            
                if (obj.orth){
                    thisOrf.attr("id", gid+'_'+co)
                        .attr('end', obj.end+pp1).attr('len', obj.L)
                        .on("click", function(){ alignOrth(co) })
                } else {
                    thisOrf.on("click", function(){ locateORF(obj.cid) })
                }

            // xpath
                if (!isUv){ return }
                var points = [];
                if (xpath[co]){ points = xpath[co] }
                points.push([start,end]);
                xpath[co] = points;
            });
            
            start0 = pp2
		});
		myStrain.attr("len", start0);

		//draw axis according to B31
		if (gid==100) {
		    var tick_mark = svg.append("g").attr("id", "tick_mark");
			drawAxis(tick_mark, r, start0, 200, 5, 1000)
		}
    });

    //path
    uvShade(xpath);

    //draw guideline
    guideline = svg.append("line").attr("class","guideline")
    			.attr("y1",2).attr("y2",hUnit*gidC.length+yEdge+5)
    			.attr("x1","-50px").attr("x2","-50px");

    $("#showSynteny").scrollLeft(0)
}

function uvShade(p){
    var pp = {};
    $.each(p, function(co,xs){
        var p1 = xs.map(function(x){ return x[0]}),
            p2 = xs.map(function(x){ return x[1]}).reverse();
        pp[co] = $.merge(p1,p2)
    });
    $.each(pp, function(co,xs){d3.select('#s'+co).attr('d', pathFun(xs))})
}


function strainTbl() {
    var title = ['Species','Geography','Source','ospC','Download','Citation'],
        posX = [44,147,247.5,316,379,444.5],
        pX = [8,92,215],
        svg_w = posX[posX.length-1] + 29.5;
    
    var svg = d3.select('#strain_tbl').append("svg").attr("width", svg_w).attr("height", bgStartL+hUnitL*gids.length);
    var svgTitle = svg.append("g").attr("class", "strain_title");

    var bg1 = svg.append("g").attr("class", "bg1"),
    	bg2 = svg.append("g").attr("class", "bg2");

    bg1.append("rect").attr("x", 0). attr("y", bgStartL).attr("width", svg_w).attr("height", hUnitL*nG1);
    bg2.append("rect").attr("x", 0). attr("y", bgStartL+hUnitL*nG1).attr("width", svg_w).attr("height", hUnitL*nG2);

    for (var j=0; j<title.length; j++){
        svgTitle.append("text").attr("x", posX[j]).attr("y", bgStartL-12).text(title[j])
    }
    
    var base = bgStartL + hUnitL;
    var text = svg.append("g").attr("class", "strainFont");
    $.each(gids, function(i,g) {
        svg.append("line").attr("x1", 0).attr("x2", svg_w).attr("y1", base).attr("y2", bgStartL+hUnitL*(i+1));

        var obj = strainData[g];
        for (var j=0; j<title.length; j++){
            var txt = text.append("text").attr("x", j>2? posX[j] : pX[j]).attr("y", base-(j<4? 4 : 1));
            
            if (!j){ txt.text(taxoData[obj.tid]).attr("class", "ita") }
            else if (j==1){ txt.text(geoData[obj.geo_id][0]).attr("class", 'geo'+obj.geo_id) }
            else if (j==2 && obj.source){
                txt.text(obj.source);
                if (/\./.test(obj.source) || obj.source=='Dermacentor'){ txt.attr("class", "ita") } }
            else if (j==3 && obj.ospc){ txt.text(obj.ospc).attr("class", "anchM") }
            else if (j==4 || j==5 && obj.citation){
                var lk = j==4? '#' :  dbLink.pubmed + obj.citation.join(',');
                var tx = txt.append("svg:a").attr("xlink:href", lk).html('&bull;').attr("class","strainBull");
                if (j==4){
                    tx.on("mouseover", function(){downFas(i,g)})
                      .attr("id", 'fa'+g)
                } else {
                    tx.on("mouseover", function(){showPub(i,obj.citation)})
                      .on("mouseout", function(){hidePub()})
                      .attr("target", "_blank")
                }
            }
        }
        base += hUnitL
    }) 
}
function showPub(i,acc) {
    var pubThing = '';
    for (var j=0; j<acc.length; j++){
        if (j){ pubThing += '<br><br>' }
        var pub =  pubData[acc[j]],
            string = pub[0] + ' et al.<br>' + pub[1] + '<br>' + pub[2] + ' ' + pub[3];
        pubThing += string
    }
    $('#pub').css("top", (i*16+55)+"px").html(pubThing).show()
}
function hidePub() { $('#pub').hide() }


var downFile = ['nuc', 'pep'];
function downFas(i,g) {
	$('#geoName').html(strainData[g].strain);
	$('#fasta').css("top", ((i>40? 40 : i)*16+70) + "px").show();
	d3.selectAll('#strain_tbl a').classed('act', false);
	d3.select('#fa'+g).classed('act', true);
    $.each($('#fasta a'), function(n, elem){
        $(elem).attr("href", "Download/"+g+'.'+downFile[n])
    })
}


function repliconTbl() {
    var wCell = 24, sp_l = 16, hContigs = 72,
        svg_w = wCell*group_id.length;
    var svgT = d3.select('#conTitle').append("svg").attr("width", svg_w).attr("height", hContigs);
    
    svgT.append("rect").attr("x", 0). attr("y", 0).attr("width", wCell*3-3).attr("height", hContigs).attr("class","grayFill");

    $.each(group_id, function(i,id) {
		svgT.append("text")
		    .attr("x", 0).attr("y", 0)
	    	.attr("transform", "translate(" + (wCell*i+15) + "," + hContigs + "),rotate(-90)")
	    	.attr("id", 'plasmid'+id)
	    	.text(groupName[id][0] + (groupName[id][1]? (' (' + groupName[id][1] + ')') : ''))
    });

    var svg = d3.select('#genome_tbl').append("svg").attr("width", svg_w).attr("height", hUnitR*gids.length);

    var bg1 = svg.append("g").attr("class", "bg1"),
    	bg2 = svg.append("g").attr("class", "bg2");

    bg1.append("rect").attr("x", 0). attr("y", 0).attr("width", svg_w).attr("height", hUnitR*nG1);
    bg2.append("rect").attr("x", 0). attr("y", hUnitR*nG1).attr("width", svg_w).attr("height", hUnitR*nG2);
    
    $.each(gids, function(i,g){
        svg.append("line").attr("x1", 0).attr("x2", svg_w).attr("y1", hUnitR*(i+1)).attr("y2", hUnitR*(i+1));
		var con = acc_genome[g];
        $.each(group_id, function(j,id){
    		if (con[id]) {
                var ar = con[id][0];
                svg.append("svg:a").attr("target", "_blank").attr("xlink:href", dbLink.ncbi + con[id].join(','))
                    .append("text")
                    .attr("x", wCell*(j+0.5))
                    .attr("y", hUnitR*(i+1)+2)
                    .attr("class", ar)
                    .html($.inArray(con[id][0], conFuse) == -1? "&#9642" : "&#9643")
                    .on("mouseover", function(){showPos(g, id, ar)})
				    .on("mouseout", function(){hidePos()})
            }
        })
    })
}
function showPos(strain, gpid, con) {
    $('#gtree #' + strain + ', #plasmid'+gpid).css("fill", "darkorange");
    d3.selectAll('.'+con).classed("darkorange", true)
}
function hidePos() {
    $('#conTitle text, #gtree text').css("fill", "black");
    d3.selectAll('#genome_tbl text').classed("darkorange", false)
}

function draw_ptree(data){
    var wLimit=tree_wPls-63, gap=6, lowerYbound=base0+4, txtOffset=4;
    
    var r = wLimit/data.w_tree;

    var svg = d3.select('#ptree').append("svg").attr("width", tree_wPls).attr("height", hUnitP*gidP.length+base0),
        line_tree = svg.append("g").style("stroke", "black");

    // draw tree:
    $.each(data.nodes, function(index,node){
		var x2 = node.xcoord*r,
			yy = node.ycoord*hUnitP + lowerYbound;

        if (node.branch_length){
            var x1 = x2-node.branch_length*r;
            line_tree.append("line").attr("x1", x1).attr("x2", x2).attr("y1", yy).attr("y2", yy)
        }

		if (node.is_Leaf) {
			svg.append("text").attr("x",wLimit+gap).attr("y",yy+txtOffset).text(strainData[node.id].strain);
		} else {
            x2 = x2? x2 : 0.5;
	    	line_tree.append("line").attr("x1", x2).attr("x2", x2).attr("y1", node.descendent_ycoord[0]*hUnitP+lowerYbound).attr("y2", node.descendent_ycoord[1]*hUnitP+lowerYbound);
		}
    })
}

function draw_ctree(data){
    var wLimit = tree_wCore - 65,
        lowerXbound=4, lowerYbound = 8,
		txtOffset=5,
        gap=6;
    
    var r = wLimit/data.w_tree;

    var svg = d3.select('#trees').append("svg")
            .attr("width", tree_wCore).attr("height", hUnit*gidC.length).style("background-color", "#dae2d6");
    
    var line_bg = svg.append("g").attr("class", "line_bg");
    for (var i=1; i<=gidC.length; i++) {
        line_bg.append("line").attr("x1", wLimit+8).attr("x2", tree_wCore).attr("y1", hUnit*i+2).attr("y2", hUnit*i+2)
    }
    
    var line_tree = svg.append("g").style("stroke", "black");

    // draw tree:
    $.each(data.nodes, function(index,node){
		var x2 = node.xcoord*r + lowerXbound,
			yy = node.ycoord*hUnit + lowerYbound;

        if (node.branch_length){
            var x1 = x2-node.branch_length*r;
            line_tree.append("line").attr("x1", x1).attr("x2", x2).attr("y1", yy).attr("y2", yy)    
        }

		if (node.is_Leaf) {
			svg.append("text").attr("x",wLimit+lowerXbound+gap).attr("y",yy+txtOffset).text(strainData[node.id].strain)
		} else {
            x2 = x2? x2 : 0.5;
	    	line_tree.append("line").attr("x1", x2).attr("x2", x2).attr("y1", node.descendent_ycoord[0]*hUnit+lowerYbound).attr("y2", node.descendent_ycoord[1]*hUnit+lowerYbound)
        }
    })
}

function draw_tree(container) {
    var svg_w, wLimit, lowerXbound=4, gap=6, sp_l, lowerYbound, bgStart, txtOffset;
    if (container=='#ltree') {
		svg_w = 550;
		wLimit = svg_w-74;
		sp_l = hUnitL;
		lowerYbound = 37;
		bgStart = bgStartL;
		txtOffset = 3
    } else {
        svg_w = 240;
        wLimit = svg_w - 63;
        sp_l = hUnitR;
        lowerYbound = 8;
		bgStart = 0;
		txtOffset = 2
    }

    var r = wLimit/treeData.w_tree;

    var svg = d3.select(container).append("svg").attr("width", svg_w).attr("height", bgStart+sp_l*gids.length),
        bg1 = svg.append("g").attr("class", "bg1"),
    	bg2 = svg.append("g").attr("class", "bg2"),
        line_bg = svg.append("g").attr("class", "line_bg"),
        line_tree = svg.append("g").style("stroke", "black");

    bg1.append("rect").attr("x", 0). attr("y", bgStart).attr("width", svg_w).attr("height", sp_l*nG1);
    bg2.append("rect").attr("x", 0). attr("y", bgStart+sp_l*nG1).attr("width", svg_w).attr("height", sp_l*nG2);

    for (var i=1; i<=gids.length; i++) {
        line_bg.append("line").attr("x1", wLimit+8).attr("x2", svg_w).attr("y1", bgStart+sp_l*i).attr("y2", bgStart+sp_l*i)
    }
    
    // draw tree:
    var y_group = [], x_group = [];
    $.each(treeData.nodes, function(index,node){
		var x2 = node.xcoord*r + lowerXbound,
			yy = node.ycoord*sp_l + lowerYbound,
            x1 = x2;

        if (node.branch_length){
            x1 -= node.branch_length*r;
            line_tree.append("line").attr("x1", x1).attr("x2", x2).attr("y1", yy).attr("y2", yy)    
        }

		if (node.parent_is_root) { y_group.push(yy); x_group.push((x1+x2)/2) }
                                                
		if (node.is_Leaf) {
			var t = svg.append("text").attr("x",wLimit+lowerXbound+gap).attr("y",yy+txtOffset).text(strainData[node.id].strain);
			if (container=='#ltree') { t.attr("class", 'geo'+strainData[node.id].geo_id) }
	    	else { t.attr("id", node.id) }
		} else {
	    	line_tree.append("line").attr("x1", x2).attr("x2", x2).attr("y1", node.descendent_ycoord[0]*sp_l+lowerYbound).attr("y2", node.descendent_ycoord[1]*sp_l+lowerYbound);
		}
    });

    if (container == '#ltree') {
        var text = svg.append("g").style("text-anchor", "middle");
	    //draw scale:
	 	var u = treeData.unit,
			scaleX = 50,
			scaleY = sp_l*gids.length-7;
		line_tree.append("line").attr("x1", scaleX).attr("x2", scaleX+u*r).attr("y1", scaleY).attr("y2", scaleY);
		line_tree.append("line").attr("x1", scaleX).attr("x2", scaleX).attr("y1", scaleY-1.5).attr("y2", scaleY+1.5);
		line_tree.append("line").attr("x1", scaleX+u*r).attr("x2", scaleX+u*r).attr("y1", scaleY-1.5).attr("y2", scaleY+1.5);

		text.append("text")
            .attr("x", scaleX+u/2*r).attr("y", scaleY+15).text(u + ' base sub/site')
            .style("font-size", "10px");

		//write group name
        text.append("text").attr("x",x_group[1]).attr("y",y_group[1]-8).text('Lyme-Disease Group');
        text.append("text").attr("x",x_group[0]).attr("y",y_group[0]-8).text('Relapsing-Fever Group')
    }
}

function drawAxis(tick_mark, ratio, l, unit, sub, k) {
    var tick = tick_mark.append("g").attr("class", "tick"),
    	mark = tick_mark.append("g").attr("class", "mark");

    tick.append("line").attr("x1", 1/ratio).attr("x2", l/ratio).attr("y1", base_tick).attr("y2", base_tick);
    tick.append("line").attr("x1", 1/ratio).attr("x2", 1/ratio).attr("y1", base_tick-hTick/2).attr("y2", base_tick+hTick).style("stroke-width", "1.6px");
    mark.append("text").attr("x",1/ratio).attr("y", base_tick-5).text(1).style("text-anchor", "start");

    for (i=1; i<=parseInt(l/unit); i++){
    	var x = i*unit/ratio;
		tick.append("line").attr("x1", x).attr("x2", x).attr("y1", i%5? base_tick : base_tick-hTick/2).attr("y2", base_tick+(i%5? hTick/2 : hTick));
		if (i%sub==0) { mark.append("text").attr("x", x).attr("y", base_tick-5).text(i*unit/k + (k==1000? 'k' : (k==1000000? 'M' : ''))) }
    }
}

function drawLeg(){
    var svg_h = 19, svg_w=214, base=8, recH=4, recW=17;

    var recCls=['tnStroke dotStroke', 'cidna'],
        lTxt = ['pseudogene or sequence error', 'RNA']

    var svg = d3.select('#legORF').append("svg").attr("width", svg_w).attr("height", svg_h);
    var rect = svg.append("g").style("fill", "none"),
        txt_leg = svg.append("g").attr("class", "txt_leg");

    var xPos=7;
    for (var i=0; i<2; i++){
		rect.append("rect")
            .attr("x", xPos)
            .attr("y", i? base+1 : base)
            .attr("width", recW)
            .attr("height", i? recH-2 : recH)
            .attr("class", recCls[i]);
 		txt_leg.append("text").attr("x", xPos+21).attr("y", base+5).text(lTxt[i]);
        xPos += recW + 143
    }
}

function drawSigLeg(){
    var svg_h = 19, svg_w=189, base=9, recH=4, recW=20;

    var recCls=['transbarSig', 'transbarNo'],
        lTxt = ['signifiicant', 'non-significant']

    var svg = d3.select('#legSig').append("svg").attr("width", svg_w).attr("height", svg_h);
    var txt_leg = svg.append("g").attr("class", "txt_leg");
    
    var xPos=7;
    for (var i=0; i<2; i++){
		svg.append("line")
            .attr("x1", xPos).attr("x2", xPos+recW)
            .attr("y1", base).attr("y2", base)
            .attr("class", recCls[i]);
 		txt_leg.append("text").attr("x", xPos+24).attr("y", base+3).text(lTxt[i]);
        xPos += recW + 73
    }
}

function numberWithCommas(x) { return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") }

function draw_treeP(nodes, w_tree, unit) {
    var w_limit = wTreeP - 4;
    var r = w_limit/w_tree;
    var lowerXbound = 1;
    var dot_size = 3;
    var gap = 4;

    var svg = d3.select('#treeP').append("svg").attr("width", lowerXbound+w_limit+dot_size).attr("height", lowerYbound+(idsP.length+2)*spH);
    var line_tree = svg.append("g").style("stroke", "black");
    var line_gap  = svg.append("g").style("stroke", "lightgray").style("stroke-dasharray", ("4,2"));

    // draw tree:
    $.each(nodes, function(index,node){
        if (node.branch_length){
            line_tree.append("line")
                .attr("x1", node.xcoord*r+lowerXbound)
                .attr("x2", (node.xcoord-node.branch_length)*r+lowerXbound)
                .attr("y1", node.ycoord*spH+lowerYbound)
                .attr("y2", node.ycoord*spH+lowerYbound)            
        }

		if (node.is_Leaf) {
	    	line_gap.append("line")
	    			.attr("x1", node.xcoord*r+lowerXbound+dot_size).attr("x2", w_limit+lowerXbound)
	    			.attr("y1", node.ycoord*spH+lowerYbound).attr("y2", node.ycoord*spH+lowerYbound);
		} else {
	    	line_tree.append("line")
	    			 .attr("x1", node.xcoord*r+lowerXbound)
	    			 .attr("x2", node.xcoord*r+lowerXbound)
	    			 .attr("y1", node.descendent_ycoord[0]*spH+lowerYbound)
	    			 .attr("y2", node.descendent_ycoord[1]*spH+lowerYbound);
		}
    });

    //draw scale:
    var x = lowerXbound+15;
    var y = lowerYbound+(idsP.length+1)*spH-6;
    line_tree.append("line").attr("x1", x).attr("x2", x+unit*r).attr("y1", y).attr("y2", y);
    line_tree.append("line").attr("x1", x).attr("x2", x).attr("y1", y-1.5).attr("y2", y+1.5);
    line_tree.append("line").attr("x1", x+unit*r).attr("x2", x+unit*r).attr("y1", y-1.5).attr("y2", y+1.5);
    
    var legPTree = svg.append("g").style("font-size", "10px");
    legPTree.append("text").attr("x", x+unit*r+gap+3).attr("y", y+3).text(unit + ' base sub/site');
    legPTree.append("text").attr("x", x).attr("y", y+18).text('(>=80% support for branches)')
}

function blast(){
    $('#newsBlast').hide();
    var seq = $('#sequence').val();
    
//sanity check
    seq = seq.replace(/ /g,'').replace(/-/g,'');
    if (!seq) { return }
    var id, ss='';
    var lines = seq.split('\n');
    $.each(lines, function(i,line){
        if (!line){ return }
        if (!id){
            if (/^>/.test(line)){
                id=line.replace(/^>/, '');
                if (!id){ id = 'query' }
            }
            else { id = 'query'; ss += line }
        } else {
            ss += line
        }
    });
    if (!ss){ return }

    if (!/^[ACDEFGHIKLMNPQRSTVWY]*$/i.test(ss)){ $('#newsBlast').html('The sequence contains non-standard codes.').show() }
    
    $.ajax({
		url: "cgi-bin/runblast.pl",
        dataType: "json",
        data: 'ev=' + $('#eValue').val() + '&gid=' + $('#strainBlast').val() + '&seq=' + ss,
		error: function(){ alert('AJAX error.') },
		success: function(data) { showBlast(data, ss.length) }
    })
}

function showBlast(data, lenQ){
    if (data.rna){ $('#newsBlast').html('This is RNA sequence.').show(); return }
    else if (data.no){ $('#blastRes').html('<h3>No hits found</h3>'); return }
    $('#blastRes').empty();

    for (var i=0; i<data.length; i++){
        var obj = data[i];
        var gapQ = (obj.seqQ[2].match(/-/g) || []).length,
            gapS = (obj.seqS[2].match(/-/g) || []).length;
        
        var hitLocus = '';
        if (!obj.orfs){ hitLocus = obj.locus }
        else {
            var id = obj.orfs.shift();
            var cont = groupName[id][0];
            if (id==31){ cont = groupName[25][0] + '|' + cont }
            
            $.each(obj.orfs, function(n,o){ hitLocus += (n? ' - ' : '') + (o[0]? o[0] : (n? 'end' : 'start')) })
        }
        
        var res = '<div>';
        res += '<p><span>Subject: <b>';
        if (!i && (!obj.orfs || obj.orfs.length==1)){ res += '<a href="#" ondblclick="getParalog(0,0,\'' + hitLocus  + '\')" onMouseOver="showInfo()" onMouseOut="hideInfo()">'}
        res += hitLocus;
        if (!i){ res += '</a>'}
        res += '</b>';
        if (cont){ res += ' at <b>' + cont + '</b>' }
        res += '</span>';

        res += '<span>E value: ' + obj.ev + '</span>';
        res += '<span>Positives: ' + Math.round(obj.Positives*100) + '%</span>';
        res += '<span>Identities: ' + Math.round(obj.Identities*100) + '%</span>';
//        res += '<span>Gaps: ' + obj.Gaps + '</span></p>';
        res += '<span>Gaps in Query: ' + gapQ + '</span>';
        res += '<span>Gaps in Subject: ' + gapS + '</span></p>';

        res += '<table><tr>';
        res += '<td class="tickSeq"><br>Query<br><br>Subject</td>';
        
        res += '<td><div class="blastSeqArea">';
        res += '<span class="tickSeq">' + obj.seqQ[0] + Array(obj.seqQ[1]-obj.seqQ[0]+gapQ-(obj.seqQ[0].length-1)).join("&nbsp;") + obj.seqQ[1] + '</span><br>';
        res += '<span class="blastSeq">' + obj.seqQ[2] + '<br>' + modiSeqH(obj.seqH) + '<br>' + obj.seqS[2] + '</span><br>';
        res += '<span class="tickSeq">' + obj.seqS[obj.rev?1:0] + Array(obj.seqS[1]-obj.seqS[0]+gapS-(obj.seqS[0].toString().length-1)).join("&nbsp;") + obj.seqS[obj.rev?0:1] + '</span><br><br>';
        res += '</div></td>';
        
        res += '<td class="blastThumb" id="bt_' + i + '"></td>';
        res += '</tr></table>';
        res += '</div>'
        $('#blastRes').append(res);

        drawBlastThumb(i, obj.lenS, obj.seqS[0], obj.seqS[1], lenQ, obj.seqQ[0], obj.seqQ[1], obj.rev, obj.locus, obj.orfs);
    }
}

function modiSeqH(seq){ return seq.replace(/[A-Y]/ig,'|').replace(/\s/g, '&nbsp;') }
function showInfo(){ $('#infoBlast').show() }
function hideInfo(){ $('#infoBlast').hide() }

var wBlast = 260;
function drawBlastThumb(i, lenS, s1, s2, lenQ, q1, q2, rev, hitLocus, orfs){
    if (orfs){
        var s0 = orfs[0]? orfs[0][1] : 0;
        s1 -= s0;
        s2 -= s0;
        lenS = (orfs[orfs.length-1]? orfs[orfs.length-1][2] : lenS) - s0
    }
    var p0=Math.max(q1,s1),
        len = p0 + Math.max(lenQ-q1, lenS-s1),
        mar = 16;
        ratioB = len/(wBlast-mar*2),
        baseQ=23, baseS=56, hRect=5, lArr=16,
        yv = [baseS-hRect, baseS-hRect, baseS-1, baseS+1, baseS+hRect, baseS+hRect];
    
    var xscale = d3.scale.linear().domain([0, len]).range([mar, wBlast-mar]),
        off = xscale(p0);

    var svg = d3.select("#bt_"+i).append("svg")
            .attr("width",wBlast).attr("height", 93);

    var unmatch = svg.append("g").attr("class", "blastLine blastUnmatch"),
        orfRect = svg.append("g").attr("class", "blastOrf"),
        match = svg.append("g").attr("class", "blastLine blastMatch"),
        guide = svg.append("g").attr("class", "blastUnmatch");
    
    if (!orfs){
        var start = off - s1/ratioB,
            end   = off + (lenS-s1)/ratioB;
        var xv = [start, end-lArr, end, end, end-lArr, start],
            points = $.map(xv, function(x,i){ return x + ',' + yv[i]}).join(' ');
        
        orfRect.append("polygon").attr("points", points);
        svg.append("text").attr("x", off+(lenS-s1*2)/ratioB/2).attr("y",72).text(hitLocus)
    } else {
        $.each(orfs, function(n,o){
            if (!o){ return }
            var start = off+(o[1]-s0-s1)/ratioB,
                end   = off+(o[2]-s0-s1)/ratioB;
            var xv = o[3]? [start, end-lArr, end, end, end-lArr, start] : [end, start+lArr, start, start, start+lArr, end],
                points = $.map(xv, function(x,i){ return x + ',' + yv[i]}).join(' ');

            orfRect.append("polygon").attr("points", points)
            svg.append("text").attr("x", start+(o[2]-o[1])/ratioB/2).attr("y",72).text(o[0])
        })
    }

    unmatch.append("line").attr("x1", off-q1/ratioB).attr("x2", off+(lenQ-q1)/ratioB).attr("y1", baseQ).attr("y2", baseQ);
    unmatch.append("line").attr("x1", off-s1/ratioB).attr("x2", off+(lenS-s1)/ratioB).attr("y1", baseS).attr("y2", baseS);
    
    guide.append("line").attr("x1", off).attr("x2", rev? (off+(s2-s1)/ratioB) : off).attr("y1", baseQ).attr("y2", baseS);
    guide.append("line").attr("x1", off+(q2-q1)/ratioB).attr("x2", rev? off : off+(s2-s1)/ratioB).attr("y1", baseQ).attr("y2", baseS);

    match.append("line").attr("x1", off).attr("x2", off+(q2-q1)/ratioB).attr("y1", baseQ).attr("y2", baseQ);
//    match.append("line").attr("x1", off).attr("x2", off+(s2-s1)/ratioB).attr("y1", baseS).attr("y2", baseS);
    
//draw axis
    var xAxis = d3.svg.axis().scale(xscale).orient("bottom").ticks(5).tickSize(-3, 0)
                .tickFormat(function(f){ return f + (s0? s0 : 0) });
    svg.append("g").attr("class", "xaxis").attr("transform", "translate(" + (p0-s1)/ratioB + ",82)").call(xAxis);
}

function initMap(){
    var lat1=0, lat2=0, lon1=0, lon2=0;
	$.each(geoData, function(i,g){
		if (!g[1] || !g[2]){ return }
		if (g[1]<lat1){ lat1 = g[1] }
		if (g[1]>lat2){ lat2 = g[1] }
		if (g[2]<lon1){ lon1 = g[2] }
		if (g[2]>lon2){ lon2 = g[2] }
	});

	var map = new google.maps.Map(document.getElementById('map'), {
		disableDefaultUI:true,
		scaleControl:true,
		center: new google.maps.LatLng((lat1+lat2)/2, (lon1+lon2)/2),
		zoom:1
    });

	$.each(geoData, function(geo_id,obj){
		col = "red";
	    var marker = new google.maps.Circle({
					center: new google.maps.LatLng(obj[1],obj[2]),
					radius: 150000,
					strokeColor: col,
					strokeWeight: 1,
					fillColor: col,
//					fillOpacity: 0.4,
					title: obj[0],
					map: map
        });

		google.maps.event.addListener(marker,'mouseover',function(){
            this.getMap().getDiv().setAttribute('title',this.get('title'))
			d3.selectAll('#ltree text, #strain_tbl text').classed('geoCol', false);
			d3.selectAll('.geo'+geo_id).classed('geoCol', true);
         });
		google.maps.event.addListener(marker,'mouseout',function(){
            this.getMap().getDiv().removeAttribute('title');
			d3.selectAll('#ltree text, #strain_tbl text').classed('geoCol', false);
         });
//		google.maps.event.addListener(marker,'click',function(){})
	})
}
