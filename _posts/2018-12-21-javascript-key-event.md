---
layout: post
title: "[JavaScript]  input 박스에 숫자만 입력 받기(자릿수 제한)"
tags: [JavaScript, JQuery, keydown, keyup, keypress, blur, change, event]
categories: [JavaScript]
subtitle: "input box + 자릿수 제한하는 직접 만든 라이브러리"
feature-img: "md/img/thumbnail/java-script-logo.png"
thumbnail: "md/img/thumbnail/java-script-logo.png"
excerpt_separator: <!--more-->
sitemap:
changefreq: daily
priority: 1.0
---

<!--more-->

# input 입력 숫자만 + 소수점 자릿수 제한 라이브러리

---

### 들어가기전

input의 값을 숫자 자릿수까지 체크하여 입력 제한하는 라이브러리를 찾던 중 마땅한 것이 없어 직접 개발하게 되었다.

`digitNumber` 전체 소스 코드는 [github](https://github.com/gmun/digit-number-javascript-key-event/blob/master/digitNumber.js)에서 확인할 수 있다.

### Version

- JQuery 3.3.1

``` javascript
<script src="https://code.jquery.com/jquery-3.3.1.js" integrity="sha256-2Kok7MbOyxpgUVvAk/HJ2jigOSYS2auK4Pfzbm7uH60=" crossorigin="anonymous"></script>
```

---

### 사용법

#### 1. 타겟 지정

먼저 적용시킬 테이블의 id 값을 지정한다.

- 반드시 id로 지정한다.

``` html
<table id="default">
    ...
</table>
```

#### 2. 라이브러리 적용

자바스크립트에서 `$(해당 table의 id).digitNumber();`로 라이브러리를 호출한다.

- `$("#default").digitNumber();`

``` javascript
$(document).ready(function(){
  $("#default").digitNumber();
});
```

#### 3.자릿수 지정

자릿수를 제한할 input 태그의 클래스 명을 `digitLimit`로 지정한다.

- `.`을 기준으로 제한할 자릿수가 결정된다.
- 앞자리가 0인 경우엔 앞자리에 0 이외의 다른 숫자를 입력을 제한한다.

```
<input type="text"  class="digitLimit1.2">
<input type="text" class="digitLimit1">
<input type="text" class="digitLimit0.2">
<input type="text" class="digitLimit.2">
<input type="text" class="digitLimit.">
<input type="text" >
```

<div>
<table id="default">
	<tr>
		<td> <label>1.2</label> </td>
		<td> <label>1</label> </td>
		<td> <label>0.2</label> </td>
		<td> <label>.2</label> </td>
    <td> <label>.</label> </td>
		<td> <label>none</label> </td>
	</tr>
	<tr>
		<td> <input type="text" class="digitLimit1.2"> </td>
		<td> <input type="text" class="digitLimit1"> </td>
		<td> <input type="text" class="digitLimit0.2"> </td>
		<td> <input type="text" class="digitLimit.2"> </td>
    <td> <input type="text" class="digitLimit."> </td>
		<td> <input type="text" > </td>
	</tr>
</table>
</div>

### 옵션) 특정 input 입력 자릿수 제한 기능 제외

테이블 안에 input 입력의 숫자와 자릿수를 제한하는 기능을 제외하기 위해선, 매개 변수에 JSON 형식의 데이터로 `notSelectors`를 전달하면 된다.

``` javascript
$(document).ready(function(){
  var options = {
    notSelectors:".notSelect, #notId" // 기능 적용을 제외할 특정 seletors
  };
  $("#notSelect").digitNumber(options);
});
```

``` html
<table id="notSelect">
    ...
    <tr>
      <td> <input type="text" > </td>
      <td> <input type="text" class="notSelect"> </td> <!-- 제외 .notSelect -->
      <td> <input type="text" id="notId"> </td> <!-- 제외  #notId -->
    </tr>
</table>
```

<div>
<table id="notSelect">
	<tr>
		<td> <label>default</label> </td>
		<td> <label>notClass</label> </td>
		<td> <label>notId</label> </td>
	</tr>
	<tr>
		<td> <input type="text" > </td>
		<td> <input type="text" class="notSelect"> </td>
		<td> <input type="text" id="notId"> </td>
	</tr>
</table>
</div>

### 옵션) 콤마 제거

#### 전체 콤마 제거

``` javascript
$(document).ready(function(){
  var options = {
    comma : false // 콤마 제거
  };
  $("#notComma").digitNumber(options); // 전체 콤마 제거
});
```

``` html
<input type="text" value="" class="digitLimit4.4">
```

<div>
<table id="notComma">
  <tr>
    <td> <label>4.4</label> </td>
  </tr>
	<tr>
		<td> <input type="text" value="" class="digitLimit4.4"> </td>
	</tr>
</table>
</div>

#### 특정 콤마 제거

특정 input에 콤마를 제거하고 싶다면 `nComma` 클래스를 추가하면된다.

``` javascript
$(document).ready(function(){
  $("#choiceNotComma").digitNumber(); // 특정 selector 콤마 제거
});
```

``` html
<input type="text" value="" class="digitLimit4.4 nComma"> <!-- nComma 추가 -->
```

<div>
<table id="choiceNotComma">
	<tr>
		<td> <label>4.4</label> </td>
		<td> <label>4.4 nComma</label> </td>
	</tr>
	<tr>
		<td> <input type="text" value="" class="digitLimit4.4"> </td>
		<td> <input type="text" value="" class="digitLimit4.4 nComma"> </td>
	</tr>
</table>
</div>

<script src="https://code.jquery.com/jquery-3.3.1.js" integrity="sha256-2Kok7MbOyxpgUVvAk/HJ2jigOSYS2auK4Pfzbm7uH60=" crossorigin="anonymous"></script>
<script type="text/javascript">
  	$(document).ready(function(){
  		$("#default").digitNumber();

  		var data = {
  					notSelectors:".notSelect, #notId"
  				};
  		$("#notSelect").digitNumber(data);

  		var data1 = {
  				comma : false
  		};
  		$("#notComma").digitNumber(data1);
  		$("#choiceNotComma").digitNumber();
  	});
</script>

<script type="text/javascript">
(function($, undefined) {
	"use strict"; // 엄격모드
	var defaults = {
			author: "Moon"
		   ,since: "2018-12-21"
		   ,project: "digitNumber"
	};

	var nk = $.digitNumber = {version: "1.0"}
	$.fn.digitNumber = function(){
		var callFn	= ""
		   ,options = {};

		for(var i in arguments){
			switch (typeof arguments[i]){
				case "string":
					callFn = arguments[i];
				break;
				case "object":
					options = arguments[i];
				break;
			}
		}

		this.each(function(i, _element) {
			var element = $(_element);
			var nKinds = new DigitNumber(element, callFn, options);
			element.data("digitNumber", nKinds);
			nKinds.render();
		});
	}

	function DigitNumber(element, callFn, options){
		var t = this;

		//export
		t.render 		= render;
		t.core			= core;
		t.initSelectors	= initSelectors(element, options);
		t.options		= options;

		function render(){
			EventManager.call(t, element);
		}

	}

	function EventManager(element){
		var t = this;

		//import
		t.core.call(t);
		t.event.call(t);

		//constract
		(function(){
			var selectors = t.initSelectors;
				setImeMode(selectors);

			t.addEvent(selectors);
		})();


		//ime-mode:disabled
		function setImeMode(selectors){
			$(selectors).css("-webkit-ime-mode", "disabled")
					    .css("-moz-ime-mode", "disabled")
					    .css("-moz-ime-mode", "disabled")
					    .css("-ms-ime-mode", "disabled")
					    .css("ime-mode", "disabled");
		}
	}

	function core(){
		var t = this;

		//import
		t.format = foramt;
		t.event  = event;
		t.digit	 = digit;
		t.regexp = regexp;
	}

	function event(){
		var t = this;

		//import
		t.core.call(t);
		t.regexp.call(t);
		t.format.call(t);
		t.digit.call(t);

		//export
		t.addEvent   = addEvent;
		t.disConnect = disConnect;

		function addEvent(selectors){
			fetchEventSource(selectors);
		}

		function disConnect(event){
			event.preventDefault();				// 현재 이벤트의 기본 동작을 중단한다.
			event.stopPropagation();			// 현재 이벤트가 상위로 전파되지 않도록 중단한다.
			event.stopImmediatePropagation();	// 현재 이벤트가 상위뿐 아니라 현재 레벨에 걸린 다른 이벤트도 동작하지 않도록 중단한다.
		}

		function fetchEventSource(selectors) {
			for(var i in selectors){
				var _selector = selectors[i];

				$(_selector).bind("change blur", function(){
					t.overLimitNumSlice(this);
					$(this).val(t.decimalComma(this, $(this).data("commaYN")));
				}).bind("keydown", function(event){
					event = event || window.event;	// chorme, ie 이벤트 구별

					//value
					var _key = event.key
					   ,_value = $(this).val();

					//위치
					var _point 	  = t.cursorPosition(this)
					   ,_dotPoint = _value.indexOf(".");

					//포함여부
					var _dotIncludeFlag = _dotPoint > -1 ? true : false;

					//자릿수 obj
					var _realLimitDigitObj = t.realLimitDigitObj(this) // 현재
					   ,_realPreDigit  = _realLimitDigitObj[0]
					   ,_realPostDigit = _realLimitDigitObj[1]
					   ,_stdLimitDigitObj  = t.stdLimitDigitObj(this)  // 기준 자릿수 obj
					   ,_stdPreDigit   = _stdLimitDigitObj[0]
					   ,_stdPostDigit  = _stdLimitDigitObj[1];

					var eventActionFlag = (_stdPreDigit === -1 && _stdPostDigit === -1) ? false : true;

					if(	   _key == "Tab"
						|| _key == "ArrowRight"
						|| _key == "ArrowLeft"
						|| _key == "Backspace"
						|| _key == "Delete"
						|| _key == "Home"
						|| _key == "End"
						){
						return;
					}else{
						if(eventActionFlag){
							if(!t.getRegexp("dotAndOnlyNumber").test(_key)){
								t.disConnect(event);
								return false;
							}

							if(_dotIncludeFlag){ //dot include
								if(_key == "."){
									t.disConnect(event);
									return false;
								}
							}else{
								if(_key == "."){
									return;
								}
							}

							if(_point > _dotPoint
								&& _dotPoint > -1){
								//dot post
								if(_stdPostDigit !== -1){
									if(_realPostDigit == _stdPostDigit
										&& !(_stdPreDigit == 0 && _key == 0)){
										t.disConnect(event);
										return false;
									}
								}
							}else{
								//dot pre
								if(_stdPreDigit !== -1){
									if(_realPreDigit == _stdPreDigit){
										if(!(_stdPreDigit == 0 && _key == 0)){
											t.disConnect(event);
											return false;
										}
									}else{
										if( _realPreDigit >_stdPreDigit){
											if(!(!_dotIncludeFlag && _key == ".")){
												t.disConnect(event);
												return false;
											}
										}
									}
								}
							}
						}
					}
					return;
				});
			}
		}


	}

	function regexp(_type){
		var t = this;
		t.getRegexp = (function(_type){
			switch (_type) {
				case "onlyNumber":	//only number
					return /^[0-9]*$/;
				break;
				case "dotAndOnlyNumber": // number or dot
					return /^[0-9]*$|\./;
				break;
				case "stdDigitClass": // stdDigitClass
					return /digitLimit/;
				break;
			}
		});
	}


	function foramt(){
		var t = this;

		//export
		t.decimalComma = decimalComma;

		function decimalComma(_selector, _comma){

			var _value = $(_selector).val().toString().replace(/[^(0-9|\.)]/gi,"") || "";

			var _realObj = _value.split(".")
			   ,_preNum  = _realObj[0] || ""
			   ,_postNum = _realObj[1] || "";

			while(/^0/gi.test(_preNum) && _preNum.length > 1){
				_preNum = _preNum.replace(/^0/gi, "");
			}

			while(/0$/gi.test(_postNum)){
				_postNum = _postNum.replace(/0$/gi, "");
			}

			var dot = _postNum === "" ? "" : ".";

			_preNum = dot === "." && _preNum == "" ? "0" : _preNum;

			if(t.options.comma !== false && !_comma){
				_preNum = comma(_preNum);
			}

			_value = _preNum + dot + _postNum;
			return _value;
		}

		function comma(_x){
			_x = decomma(_x);
			return _x.replace(/\B(?=(\d{3})+(?!\d))/g, ",");
		}

		function decomma(_x){
			return _x.toString().replace(/[^(0-9)]/gi,"");
		}
	}

	function digit(){
		var t = this;

		//import
		t.regexp.call(t);

		//export
		t.cursorPosition 	= cursorPosition;		// 현재 자릿수 반환
		t.overLimitNumSlice = overLimitNumSlice;	// 자릿수 넘어가면 자르기
		t.realLimitDigitObj = realLimitDigitObj;	// 현재 자릿수 obj
		t.stdLimitDigitObj	= stdLimitDigitObj;		// 기준 자릿수 obj

		function cursorPosition(_selector) {
			_selector = $(_selector);
		    var start = _selector[0].selectionStart
		       ,end   = _selector[0].selectionEnd
		       ,diff  = end - start;

			/*
		    if (start >= 0 && start == end) {
		        // do cursor position actions, example:
		        //'Cursor Position: ' + start
		    } else if (start >= 0) {
		        // do ranged select actions, example:
		        //'Cursor Position: ' + start + ' to ' + end + ' (' + diff + ' selected chars)'
		    }
		    */
	  	  return start;
		}

		function overLimitNumSlice(_selector){
			var _realLimitDigitObj = t.realLimitDigitObj(_selector)
			   ,_stdLimitDigitObj  = t.stdLimitDigitObj(_selector);

			var eventActionFlag = (_stdLimitDigitObj[0] === -1 && _stdLimitDigitObj[1] === -1) ? false : true;

			if(eventActionFlag){
				var _value = $(_selector).val().toString().replace(/[^(0-9)|\.]/gi,"");
				if(_value === ""){
					return;
				}
				var _stdPreDigit  = _stdLimitDigitObj[0]
				   ,_stdPostDigit = _stdLimitDigitObj[1]
				   ,realPreDigit  = _realLimitDigitObj[0]
				   ,realPostDigit = _realLimitDigitObj[1];

				var tmpObj = _value.split(".");
				var _preNumb  = tmpObj[0] || "0"
				   ,_postNumb = tmpObj[1] || "0";

				if(_stdPreDigit !== -1){
					_preNumb = _preNumb.substring(0, _stdPreDigit) || 0;
				}

				if(_stdPostDigit !== -1){
					_postNumb = _postNumb.substring(0, _stdPostDigit) || 0;
				}

				_value = _preNumb + "." + _postNumb;

				$(_selector).val(_value);
			}
		}

		function realLimitDigitObj(_selector){
			var _limitObj = [0, 0];

			var _value = $(_selector).val().toString().replace(/[^(0-9)|\.]/gi,"");

			var _tmpObj = _value.split(".");
			if(_tmpObj.length > 0){
				var _preNumb  = _tmpObj[0] || ""
				   ,_postNumb = _tmpObj[1] || "";

				//자릿수
				_preNumb  = _preNumb.length;  // 정수
				_postNumb = _postNumb.length; // 소수점

				_limitObj = [_preNumb, _postNumb];
			}

			return _limitObj;

		}

		function stdLimitDigitObj (_selector){
			var _limitObj = [-1, -1];

			var _classList = $(_selector)[0].classList;
			var _limitDigit = "";

			var _classPattern = t.getRegexp("stdDigitClass");

			for(var i=0, objLength = _classList.length; i<objLength; i++){
				var _className = _classList[i];

				if(_classPattern.test(_className)){
					_limitDigit = _className.replace(_classPattern, "");
				}
			}

			if(_limitDigit !== ""){
				_limitObj = _limitDigit.split(".");
				var _preNumb  = _limitObj[0] || -1
				   ,_postNumb = _limitObj[1] || -1;

				_limitObj = [_preNumb, _postNumb];
			}

			return _limitObj;
		}
	}

	function initSelectors(element, options){
		var _reObj = new Array();

		if(element !== undefined){
			if(element[0].id == ""){
				console.error("ERROR :: Should must you have to ID, not a class name. \nex) $('#id').digitNumber(..);");
				return _reObj;
			}

			var inputSeletors = $("#" + element[0].id + " input");

			var notSeletors = options.notSelectors || "";
				notSeletors = notSeletors.replace(/ /gi, "");

			var _notObj = notSeletors.split(",");

			for(var i=0, tot=inputSeletors.length; i<tot; i++){
				var include = true;

				$(inputSeletors[i]).data("commaYN", $(inputSeletors[i]).hasClass("nComma"));

				for(var idx in _notObj){
					var selector = _notObj[idx].substring(0, 1);
					if(selector !== ""){
						var _nS 	 = _notObj[idx].replace(new RegExp("\\" + selector), "");
						var _nSRegex = new RegExp("\\b(" + _nS + ")\\b", "g");

						var selectorStr;
						switch (selector) {
							case "#":	// id
								selectorStr = inputSeletors[i].id;
							break;

							case ".":	// class
								selectorStr = inputSeletors[i].classList.toString();
							break;
						}

						if(_nSRegex.test(selectorStr)){
							include = false;
							break;
						}
					}
				}

				if(include){
					_reObj.push(inputSeletors[i]);
				}

			}
		}
		return _reObj;
	}

})(jQuery);

</script>
