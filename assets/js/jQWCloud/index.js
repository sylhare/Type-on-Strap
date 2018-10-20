function createWordClud(){
	
	return (function(){

		var wordData = [];
		var tagList = $(".tag-list>a>p");
		
		for(var i=0; i<tagList.length; i++){
			var name 	= tagList[i].textContent || tagList[i].outerText || tagList[i].innerText;
			var weight 	= name.length;
				weight 	= weight*100 - weight*3;	// size 조절
			
			var obj = new Array();
				obj.word 	= name;
				obj.weight  = weight;
			//{word: 'Prashant', weight: 40}, {word: 'Ajinkya', weight: 11, color: 'green'}
			wordData.push(obj);
		}		

		$("#wordCloud").jQWCloud({
			
			 words		: wordData
		  //,cloud_color: 'yellow' 	// 전체 적인 font color	
			,minFont	: 40		// 최소 font size
			,maxFont	: 80		// 최대 font size
		  //,fontOffset	: 5
			,cloud_font_family	: "Source Sans Pro"	// font style
		  //,verticalEnabled	: false		// 수직정렬
			,padding_left		: 1			// 왼쪽 패딩 size
		  //,showSpaceDIV		: true		// 테두리
		  //,spaceDIVColor		: 'white'	// 테두리 색상
			,word_common_classes: 'WordClass'	// class 명
			,word_mouseEnter :function(){
				$(this).css("text-decoration","underline");
			}
			,word_mouseOut :function(){
				$(this).css("text-decoration","none");	
			}
			,word_click: function(){ 			
				var tag = $(this).text().trim();
				$("#" + tag).animate();
				window.location.href = "/tags/#" + tag;
			}		              
			,beforeCloudRender: function(){
			    date1 = new Date();
		 	}
		 	,afterCloudRender: function(){
		 		var date2 = new Date();
				//console.log("Cloud Completed in "+(date2.getTime()-date1.getTime()) +" milliseconds");
			}
		});
	})();
};
