{% if site.data.social.rss or site.theme_settings.rss %}
<li>
    <a feed.xml href="{{ site.data.social.feed.path | default: 'feed.xml' | relative_url }}"
       title="{{ site.data.language.str_rss_follow | default: 'Follow RSS feed' }}">
        <span class="fa-stack fa-lg">
            <i class="fas fa-circle fa-stack-2x"></i>
            <i class="fas fa-rss fa-stack-1x fa-inverse"></i>
        </span>
    </a>
</li>
{% endif %}

{% if site.data.social.email_address or site.theme_settings.email_address %}
<li>
    <a href="mailto:{{ site.data.social.email_address }}"
       title="{{ site.data.language.str_email }}">
		<span class="fa-stack fa-lg">
            <i class="fas fa-circle fa-stack-2x"></i>
            <i class="fas fa-envelope fa-stack-1x fa-inverse"></i>
        </span>
    </a>
</li>
{% endif %}

{% assign bh_url = "https://www.behance.net/" | append: site.data.social.behance %}
{% assign bh_on = site.data.social.behance or site.theme_settings.behance %}
{% include social/icon_partial.html isDisplayed=bh_on url=bh_url social="Behance" %}

{% assign bb_url = "https://bitbucket.org/" | append: site.data.social.bitbucket %}
{% assign bb_on = site.data.social.bitbucket or site.theme_settings.bitbucket %}
{% include social/icon_partial.html isDisplayed=bb_on url=bb_url social="Bitbucket" %}

{% assign db_url = "https://www.dribbble.com/" | append: site.data.social.dribbble %}
{% assign db_on = site.data.social.dribbble or site.theme_settings.dribbble %}
{% include social/icon_partial.html isDisplayed=db_on url=db_url social="Dribbble" %}

{% assign fb_url = "https://www.facebook.com/" | append: site.data.social.facebook %}
{% assign fb_on = site.data.social.facebook or site.theme_settings.facebook %}
{% include social/icon_partial.html isDisplayed=fb_on url=fb_url social="Facebook" %}

{% assign fk_url = "https://flickr.com/photos/" | append: site.data.social.flickr %}
{% assign fk_on = site.data.social.flickr or site.theme_settings.flickr %}
{% include social/icon_partial.html isDisplayed=fk_on url=fk_url social="Flickr" %}

{% assign gh_url = "https://github.com/" | append: site.data.social.github %}
{% assign gh_on = site.data.social.github or site.theme_settings.github %}
{% include social/icon_partial.html isDisplayed=gh_on url=gh_url social="GitHub" %}

{% assign ig_url = "https://instagram.com/" | append: site.data.social.instagram %}
{% assign ig_on = site.data.social.instagram or site.theme_settings.instagram %}
{% include social/icon_partial.html isDisplayed=ig_on url=ig_url social="instagram" %}

{% assign kb_url = "https://keybase.io/" | append: site.data.social.keybase %}
{% assign kb_on = site.data.social.keybase or site.theme_settings.keybase %}
{% include social/icon_partial.html isDisplayed=kb_on url=kb_url social="Keybase" %}

{% assign se_url = site.data.social.linkedin %}
{% assign se_on = site.data.social.linkedin or site.theme_settings.linkedin %}
{% include social/icon_partial.html isDisplayed=se_on url=se_url social="Linkedin" %}

{% assign pi_url = "https://www.pinterest.com/" | append: site.data.social.pinterest %}
{% assign pi_on = site.data.social.pinterest or site.theme_settings.pinterest %}
{% include social/icon_partial.html isDisplayed=pi_on url=pi_url social="Pinterest" %}

{% assign rd_url = "https://www.reddit.com/user/" | append: site.data.social.reddit %}
{% assign rd_on = site.data.social.reddit or site.theme_settings.reddit %}
{% include social/icon_partial.html isDisplayed=rd_on url=rd_url social="Reddit" %}

{% assign sc_url = "https://soundcloud.com/" | append: site.data.social.soundcloud %}
{% assign sc_on = site.data.social.soundcloud or site.theme_settings.soundcloud %}
{% include social/icon_partial.html isDisplayed=sc_on url=sc_url social="SoundCloud" %}

{% assign se_url = site.data.social.stack_exchange %}
{% assign se_on = site.data.social.stack_exchange or site.theme_settings.stack_exchange %}
{% include social/icon_partial.html isDisplayed=se_on url=se_url social="Stack-Exchange" %}

{% assign so_url = site.data.social.stack_overflow %}
{% assign so_on = site.data.social.stack_overflow or site.theme_settings.stack_overflow %}
{% include social/icon_partial.html isDisplayed=so_on url=so_url social="Stack-Overflow" %}

{% assign st_url = "http://steamcommunity.com/id/" | append: site.data.social.steam %}
{% assign st_on = site.data.social.steam or site.theme_settings.steam %}
{% include social/icon_partial.html isDisplayed=st_on url=st_url social="Steam" %}

{% assign tb_url = "https://" | append: site.data.social.tumblr | append: ".tumblr.com/" %}
{% assign tb_on = site.data.social.tumblr or site.theme_settings.tumblr %}
{% include social/icon_partial.html isDisplayed=tb_on url=tb_url social="Tumblr" %}

{% assign gl_url = "https://gitlab.com/" | append: site.data.social.gitlab %}
{% assign gl_on = site.data.social.gitlab or site.theme_settings.gitlab %}
{% include social/icon_partial.html isDisplayed=gl_on url=gl_url social="Gitlab" %}

{% assign tw_url = "https://twitter.com/" | append: site.data.social.twitter %}
{% assign tw_on = site.data.social.twitter or site.theme_settings.twitter %}
{% include social/icon_partial.html isDisplayed=tw_on url=tw_url social="Twitter" %}

{% assign v_url = "https://vimeo.com/" | append: site.data.social.vimeo %}
{% assign v_on = site.data.social.vimeo or site.theme_settings.vimeo %}
{% include social/icon_partial.html isDisplayed=v_on url=v_url social="Vimeo" %}

{% assign wp_url = "https://" | append: site.data.social.wordpress | append: ".wordpress.com/" %}
{% assign wp_on = site.data.social.wordpress or site.theme_settings.wordpress %}
{% include social/icon_partial.html isDisplayed=wp_on url=wp_url social="WordPress" %}

{% assign yt_url = "https://www.youtube.com/channel/" | append: site.data.social.youtube %}
{% assign yt_on = site.data.social.youtube or site.theme_settings.youtube %}
{% include social/icon_partial.html isDisplayed=yt_on url=yt_url social="Youtube" %}
