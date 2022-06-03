---
layout: post
title: Grafana, nginx reverse-proxy and Docker
description: Deploy Grafana using Docker and nginx as reverse-proxy providing health-check and more advanced features such as blocking requests.
author-id: "galera"
categories: [docker]
tags: [grafana,docker,devops,nginx]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/grafana/featured.jpg"
thumbnail: "assets/img/posts/grafana/featured.jpg"
image: "assets/img/posts/grafana/featured.jpg"
---
<p>We want to start monitoring our AWS resources using Cloudwatch. The API is awesome, however its visualisation tool sucks, that's why we choose Grafana to present the data.</p>
<p><!--more--></p>
<p>Currently, we're using ECS to deploy our services and Grafana is available as a Docker image therefore this tool is a natural fit! We checked the documentation and everything looks fine: <a href="https://grafana.com/docs/installation/docker/">https://grafana.com/docs/installation/docker</a></p>
<p>The problem comes when we tried to deploy Grafana docker image using our governance tool. This tool expects a <em>/health-check </em>endpoint to detect the status of our applications. Here we have two approaches:</p>
<ul>
<li>Try to modify the Grafana code to add this logic inside</li>
<li>Add something to the Docker container to answer this endpoint.</li>
</ul>
<p>For obvious reasons we chose the second version. Initially we were thinking on some kind of script. However we realised that the requests to Grafana need to be proxied. This is even mentioned in their documentation! <a href="https://grafana.com/docs/installation/behind_proxy/">https://grafana.com/docs/installation/behind_proxy/</a> .</p>
<h5>The solution</h5>
<p>Now that we have discarded the scripting, we chose nginx to implement the reverse proxy and we delegate the health-check endpoint to a static content under webroot of nginx.</p>
<p>Hereunder is the Dockerfile, which is self-explanatory:</p>

```docker
FROM grafana/grafana
EXPOSE 8080 8080
COPY health-check /health-check
COPY start-nginx-grafana.sh /start-nginx-grafana.sh
USER root
RUN apt-get update &amp;&amp; apt-get install -y nginx
RUN chown -R grafana:grafana /etc/nginx/nginx.conf /var/log/nginx /var/lib/nginx /start-nginx-grafana.sh
RUN chmod +x /start-nginx-grafana.sh
USER grafana
RUN cp /health-check/nginx.conf /etc/nginx/nginx.conf
ENTRYPOINT [ "/start-nginx-grafana.sh" ]
```

<p>The tricky part we found was that installing nginx required sudo permissions, however this could be easily achieved changing to the user root in the Dockerfile. Grafana service runs as grafana user, so, some permissions of files and folders of the nginx services need to be changed to the grafana user.</p>
<p>The following snippet shows the nginx.conf file:</p>

```
worker_processes  1;
pid /var/lib/nginx/nginx.pid;

events {
    worker_connections  1024;
}

http {
    include       mime.types;
    default_type  application/octet-stream;
    sendfile        on;
    keepalive_timeout  65;
    server {
        listen       8080;
        server_name  localhost;

        location /health-check {
            default_type  "application/json";
            root   /health-check;
            index  health-check.json;
        }

        location / {
            proxy_pass http://localhost:3000/;
        }

        #location /api {
        #    return 403;
        #}
    }
}
```

<p>This configuration enables the health-check endpoint to be compatible with our governance tool. Precisely nginx returns the file health-check.json in response to this endpoint. nginx proxies any other request to the Grafana instance running inside the container. The presence of the nginx reverse proxy enables the user to implement more features, like blocking the Grafana API, etc...</p>
<p>You can check the code to provide Grafana,nginx and Docker here: <a href="https://github.com/adriangalera/nginx-grafana-docker">https://github.com/adriangalera/nginx-grafana-docker</a></p>
