#include<string.h>
#include<netdb.h>
#include<arpa/inet.h>

char * h2ip(const char *hostname) {
	struct hostent *he;
	struct in_addr **addr_list;

	if( (he = gethostbyname(hostname)) == NULL) {
		herror("gethostbyname");
		return "";
	}

	addr_list = (struct in_addr **) he->h_addr_list;
	for(int i=0; addr_list[i] != NULL; i++) {
		return strdup(inet_ntoa(*addr_list[i]));
	}
	return "";
}
