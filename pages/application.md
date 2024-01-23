---
layout: page
title: Application
feature-img: "assets/img/header/tab_back.png"
position: 4
---

### Application to enroll the course

<p align="justify">
Applicants should fill the form below in English <u>prior to March 31, 2024 (AoE).</u></p><br>

<i>Up to 25 attendees will be accepted. Applicants will be selected based on the potential value of the course for them, with their motivation letter and CV serving as criteria. Selected applicants will be notified no later than April 15.</i><br><br>

<p align="justify"><u>Instructions for attachments:</u></p>

* All documents must be uploaded in PDF format
* CV (max. 2 pages): include experience, relevant publications, talks in conferences, etc.
* Motivation letter (max. 1 page): include experience with computational chemistry and/or machine learning, what the applicant hopes to gain from the course, etc.
* GEQC membership: proof of GEQC membership. For example, a screenshot of your suscription dates from the <a href='https://rseq.playoffinformatica.com/FormLogin.php'>"Member area" of the RSEQ webpage.</a>


<br><br>

<html>
  <head>
    <title>CAMLC24 registration form</title>
    <style>
      html, body {
      min-height: 100%;
      }
      body, div, form, input, select, p { 
      padding: 0;
      margin: 0;
      outline: none;
      }
      body {
      background-size: cover;
      }
      h1, h2 {
      text-transform: uppercase;
      font-weight: 400;
      }
      h2 {
      margin: 0 0 0 8px;
      }
      .main-block {
      display: flex;
      flex-direction: column;
      justify-content: center;
      align-items: center;
      height: 100%;
      padding: 25px;
      background: #eaeaea;
      }
      .left-part, form {
      padding: 25px;
      }
      .left-part {
      text-align: center;
      }
      .fa-graduation-cap {
      font-size: 72px;
      }
      form {
      background: #ffffff; 
      border: 1px solid #ccc;
      width: 50%;
      }
      .title {
      display: flex;
      align-items: center;
      margin-bottom: 20px;
      }
      .info {
      display: flex;
      flex-direction: column;
      }
      input, select {
      padding: 5px;
      margin-bottom: 30px;
      background: transparent;
      border: none;
      border-bottom: 1px solid #ccc;
      }
      select option {
      /* margin: 40px; */
      background: #ffffff;
      color: #000000;
      /* text-shadow: 0 1px 0 rgba(0, 0, 0, 0.4); */
      }
      input, textarea {
      padding: 5px;
      margin-bottom: 30px;
      background: transparent;
      border: 1px solid #ccc;
      }
      textarea::placeholder {
      color: #ccc;
      }      
      input::placeholder {
      color: #ccc;
      }
      option:focus {
      border: none;
      }
      option {
      background: black; 
      border: none;
      }
      .btn-item, button {
      padding: 10px 5px;
      margin-top: 20px;
      border-radius: 5px; 
      border: none;
      background: #26a9e0; 
      text-decoration: none;
      font-size: 15px;
      font-weight: 400;
      color: #fff;
      }
      .btn-item {
      display: inline-block;
      margin: 20px 5px 0;
      }
      button {
      width: 100%;
      }
      button:hover, .btn-item:hover {
      background: #85d6de;
      }
      .main-block {
      flex-direction: row;
      height: calc(100% - 50px);
      }
      
      }
    </style>
  </head>
  <body>
    <center>
    <form action="https://api.web3forms.com/submit" method="POST">
    <div class="title">
        <i class="fas fa-pencil-alt"></i> 
        <h2>Register here</h2>
    </div>
    <div class="info">
        <input type="hidden" name="access_key" value="47eccb18-9823-408a-b108-e0b03e4736b5">
        <!-- Personal info -->
        <input class="fname" type="text" name="name" placeholder="Full name" required>
        <input type="email" name="email" placeholder="Email" required>
        <input type="text" name="institution" placeholder="Institution" required>
        <input type="text" name="group" placeholder="Research group" required>
        <!-- Yes or No questions -->
        <p align="justify">Are you a GEQC member?</p>
        <select name="GEQC member" required>
        <option value=""></option>
        <option value="No">No</option>
        <option value="Yes">Yes</option>
        </select>
        <p align="justify">Do you have a laptop you can bring to the course?</p>
        <select name="Laptop" required>
        <option value=""></option>
        <option value="No">No</option>
        <option value="Yes">Yes</option>
        </select>
        <p align="justify">Do you have proficient understanding of English?</p>
        <select name="English level" required>
        <option value=""></option>
        <option value="No">No</option>
        <option value="Yes">Yes</option>
        </select>
        <!-- CV box -->
        <!-- <p align="justify">
        <label for="counter-input" class="label">Characters: <span id="counter-display" class="tag is-success">0/2000</span></label></p>
        <textarea name="CV" id='counter-input' placeholder="CV summary (experience, relevant publications, talks in conferences, etc.)" rows="5" cols="50" maxlength="2000" required></textarea>
        <script>
        (() => {
          const counter = (() => {
            const input = document.getElementById('counter-input'),
              display = document.getElementById('counter-display'),
              changeEvent = (evt) => {
                const charCount = evt.target.value.length;
                display.innerHTML = `${charCount}/2000`;
              },
              getInput = () => input.value,
              countEvent = () => input.addEventListener('keyup', changeEvent),
              init = () => countEvent();
            return {
              init: init
            }
          })();
          counter.init();
        })();
        </script> -->
        <!-- Motivation letter box -->
        <!-- <p align="justify">
        <label for="counter-input2" class="label">Characters: <span id="counter-display2" class="tag is-success">0/2000</span></label></p>
        <textarea name="Motivation letter" id='counter-input2' placeholder="Motivation letter (experience with computational chemistry, what the applicant hopes to gain from the course, etc.)" rows="5" cols="50" maxlength="2000" required></textarea>
        <script>
        (() => {
          const counter = (() => {
            const input = document.getElementById('counter-input2'),
              display = document.getElementById('counter-display2'),
              changeEvent = (evt) => {
                const charCount = evt.target.value.length;
                display.innerHTML = `${charCount}/2000`;
              },
              getInput = () => input.value,
              countEvent = () => input.addEventListener('keyup', changeEvent),
              init = () => countEvent();
            return {
              init: init
            }
          })();
          counter.init();
        })();
        </script> -->
        <!-- ATTACHMENTS -->
        <p align="justify">&nbsp;&nbsp;CV of the applicant (max. 2 pages)</p>
        <input type="file" name="CV" accept="application/pdf" required>
        <p align="justify">&nbsp;&nbsp;Motivation letter (max. 1 page)</p>
        <input type="file" name="Letter" accept="application/pdf" required>
        <p align="justify">&nbsp;&nbsp;GEQC membership</p>
        <input type="file" name="GEQC" accept="application/pdf">
    </div>
    <!-- <div class="checkbox">
        <input type="checkbox" name="checkbox"><span>I agree to the <a href="https://www.w3docs.com/privacy-policy">Privacy Poalicy for W3Docs.</a></span>
    </div> -->
    <button type="submit" href="/">Submit</button>
    </form></center>
  </body>
</html>