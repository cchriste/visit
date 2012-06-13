/*
 * SSH backend.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <limits.h>
#include <signal.h>

#include "putty.h"
#include "tree234.h"
#include "ssh.h"
#ifndef NO_GSSAPI
#include "sshgssc.h"
#include "sshgss.h"
#endif

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

#define SSH1_MSG_DISCONNECT                       1	/* 0x1 */
#define SSH1_SMSG_PUBLIC_KEY                      2	/* 0x2 */
#define SSH1_CMSG_SESSION_KEY                     3	/* 0x3 */
#define SSH1_CMSG_USER                            4	/* 0x4 */
#define SSH1_CMSG_AUTH_RSA                        6	/* 0x6 */
#define SSH1_SMSG_AUTH_RSA_CHALLENGE              7	/* 0x7 */
#define SSH1_CMSG_AUTH_RSA_RESPONSE               8	/* 0x8 */
#define SSH1_CMSG_AUTH_PASSWORD                   9	/* 0x9 */
#define SSH1_CMSG_REQUEST_PTY                     10	/* 0xa */
#define SSH1_CMSG_WINDOW_SIZE                     11	/* 0xb */
#define SSH1_CMSG_EXEC_SHELL                      12	/* 0xc */
#define SSH1_CMSG_EXEC_CMD                        13	/* 0xd */
#define SSH1_SMSG_SUCCESS                         14	/* 0xe */
#define SSH1_SMSG_FAILURE                         15	/* 0xf */
#define SSH1_CMSG_STDIN_DATA                      16	/* 0x10 */
#define SSH1_SMSG_STDOUT_DATA                     17	/* 0x11 */
#define SSH1_SMSG_STDERR_DATA                     18	/* 0x12 */
#define SSH1_CMSG_EOF                             19	/* 0x13 */
#define SSH1_SMSG_EXIT_STATUS                     20	/* 0x14 */
#define SSH1_MSG_CHANNEL_OPEN_CONFIRMATION        21	/* 0x15 */
#define SSH1_MSG_CHANNEL_OPEN_FAILURE             22	/* 0x16 */
#define SSH1_MSG_CHANNEL_DATA                     23	/* 0x17 */
#define SSH1_MSG_CHANNEL_CLOSE                    24	/* 0x18 */
#define SSH1_MSG_CHANNEL_CLOSE_CONFIRMATION       25	/* 0x19 */
#define SSH1_SMSG_X11_OPEN                        27	/* 0x1b */
#define SSH1_CMSG_PORT_FORWARD_REQUEST            28	/* 0x1c */
#define SSH1_MSG_PORT_OPEN                        29	/* 0x1d */
#define SSH1_CMSG_AGENT_REQUEST_FORWARDING        30	/* 0x1e */
#define SSH1_SMSG_AGENT_OPEN                      31	/* 0x1f */
#define SSH1_MSG_IGNORE                           32	/* 0x20 */
#define SSH1_CMSG_EXIT_CONFIRMATION               33	/* 0x21 */
#define SSH1_CMSG_X11_REQUEST_FORWARDING          34	/* 0x22 */
#define SSH1_CMSG_AUTH_RHOSTS_RSA                 35	/* 0x23 */
#define SSH1_MSG_DEBUG                            36	/* 0x24 */
#define SSH1_CMSG_REQUEST_COMPRESSION             37	/* 0x25 */
#define SSH1_CMSG_AUTH_TIS                        39	/* 0x27 */
#define SSH1_SMSG_AUTH_TIS_CHALLENGE              40	/* 0x28 */
#define SSH1_CMSG_AUTH_TIS_RESPONSE               41	/* 0x29 */
#define SSH1_CMSG_AUTH_CCARD                      70	/* 0x46 */
#define SSH1_SMSG_AUTH_CCARD_CHALLENGE            71	/* 0x47 */
#define SSH1_CMSG_AUTH_CCARD_RESPONSE             72	/* 0x48 */

#define SSH1_AUTH_RHOSTS                          1	/* 0x1 */
#define SSH1_AUTH_RSA                             2	/* 0x2 */
#define SSH1_AUTH_PASSWORD                        3	/* 0x3 */
#define SSH1_AUTH_RHOSTS_RSA                      4	/* 0x4 */
#define SSH1_AUTH_TIS                             5	/* 0x5 */
#define SSH1_AUTH_CCARD                           16	/* 0x10 */

#define SSH1_PROTOFLAG_SCREEN_NUMBER              1	/* 0x1 */
/* Mask for protoflags we will echo back to server if seen */
#define SSH1_PROTOFLAGS_SUPPORTED                 0	/* 0x1 */

#define SSH2_MSG_DISCONNECT                       1	/* 0x1 */
#define SSH2_MSG_IGNORE                           2	/* 0x2 */
#define SSH2_MSG_UNIMPLEMENTED                    3	/* 0x3 */
#define SSH2_MSG_DEBUG                            4	/* 0x4 */
#define SSH2_MSG_SERVICE_REQUEST                  5	/* 0x5 */
#define SSH2_MSG_SERVICE_ACCEPT                   6	/* 0x6 */
#define SSH2_MSG_KEXINIT                          20	/* 0x14 */
#define SSH2_MSG_NEWKEYS                          21	/* 0x15 */
#define SSH2_MSG_KEXDH_INIT                       30	/* 0x1e */
#define SSH2_MSG_KEXDH_REPLY                      31	/* 0x1f */
#define SSH2_MSG_KEX_DH_GEX_REQUEST               30	/* 0x1e */
#define SSH2_MSG_KEX_DH_GEX_GROUP                 31	/* 0x1f */
#define SSH2_MSG_KEX_DH_GEX_INIT                  32	/* 0x20 */
#define SSH2_MSG_KEX_DH_GEX_REPLY                 33	/* 0x21 */
#define SSH2_MSG_KEXRSA_PUBKEY                    30    /* 0x1e */
#define SSH2_MSG_KEXRSA_SECRET                    31    /* 0x1f */
#define SSH2_MSG_KEXRSA_DONE                      32    /* 0x20 */
#define SSH2_MSG_USERAUTH_REQUEST                 50	/* 0x32 */
#define SSH2_MSG_USERAUTH_FAILURE                 51	/* 0x33 */
#define SSH2_MSG_USERAUTH_SUCCESS                 52	/* 0x34 */
#define SSH2_MSG_USERAUTH_BANNER                  53	/* 0x35 */
#define SSH2_MSG_USERAUTH_PK_OK                   60	/* 0x3c */
#define SSH2_MSG_USERAUTH_PASSWD_CHANGEREQ        60	/* 0x3c */
#define SSH2_MSG_USERAUTH_INFO_REQUEST            60	/* 0x3c */
#define SSH2_MSG_USERAUTH_INFO_RESPONSE           61	/* 0x3d */
#define SSH2_MSG_GLOBAL_REQUEST                   80	/* 0x50 */
#define SSH2_MSG_REQUEST_SUCCESS                  81	/* 0x51 */
#define SSH2_MSG_REQUEST_FAILURE                  82	/* 0x52 */
#define SSH2_MSG_CHANNEL_OPEN                     90	/* 0x5a */
#define SSH2_MSG_CHANNEL_OPEN_CONFIRMATION        91	/* 0x5b */
#define SSH2_MSG_CHANNEL_OPEN_FAILURE             92	/* 0x5c */
#define SSH2_MSG_CHANNEL_WINDOW_ADJUST            93	/* 0x5d */
#define SSH2_MSG_CHANNEL_DATA                     94	/* 0x5e */
#define SSH2_MSG_CHANNEL_EXTENDED_DATA            95	/* 0x5f */
#define SSH2_MSG_CHANNEL_EOF                      96	/* 0x60 */
#define SSH2_MSG_CHANNEL_CLOSE                    97	/* 0x61 */
#define SSH2_MSG_CHANNEL_REQUEST                  98	/* 0x62 */
#define SSH2_MSG_CHANNEL_SUCCESS                  99	/* 0x63 */
#define SSH2_MSG_CHANNEL_FAILURE                  100	/* 0x64 */
#define SSH2_MSG_USERAUTH_GSSAPI_RESPONSE               60
#define SSH2_MSG_USERAUTH_GSSAPI_TOKEN                  61
#define SSH2_MSG_USERAUTH_GSSAPI_EXCHANGE_COMPLETE      63
#define SSH2_MSG_USERAUTH_GSSAPI_ERROR                  64
#define SSH2_MSG_USERAUTH_GSSAPI_ERRTOK                 65
#define SSH2_MSG_USERAUTH_GSSAPI_MIC                    66

/*
 * Packet type contexts, so that ssh2_pkt_type can correctly decode
 * the ambiguous type numbers back into the correct type strings.
 */
typedef enum {
    SSH2_PKTCTX_NOKEX,
    SSH2_PKTCTX_DHGROUP,
    SSH2_PKTCTX_DHGEX,
    SSH2_PKTCTX_RSAKEX
} Pkt_KCtx;
typedef enum {
    SSH2_PKTCTX_NOAUTH,
    SSH2_PKTCTX_PUBLICKEY,
    SSH2_PKTCTX_PASSWORD,
    SSH2_PKTCTX_GSSAPI,
    SSH2_PKTCTX_KBDINTER
} Pkt_ACtx;

#define SSH2_DISCONNECT_HOST_NOT_ALLOWED_TO_CONNECT 1	/* 0x1 */
#define SSH2_DISCONNECT_PROTOCOL_ERROR            2	/* 0x2 */
#define SSH2_DISCONNECT_KEY_EXCHANGE_FAILED       3	/* 0x3 */
#define SSH2_DISCONNECT_HOST_AUTHENTICATION_FAILED 4	/* 0x4 */
#define SSH2_DISCONNECT_MAC_ERROR                 5	/* 0x5 */
#define SSH2_DISCONNECT_COMPRESSION_ERROR         6	/* 0x6 */
#define SSH2_DISCONNECT_SERVICE_NOT_AVAILABLE     7	/* 0x7 */
#define SSH2_DISCONNECT_PROTOCOL_VERSION_NOT_SUPPORTED 8	/* 0x8 */
#define SSH2_DISCONNECT_HOST_KEY_NOT_VERIFIABLE   9	/* 0x9 */
#define SSH2_DISCONNECT_CONNECTION_LOST           10	/* 0xa */
#define SSH2_DISCONNECT_BY_APPLICATION            11	/* 0xb */
#define SSH2_DISCONNECT_TOO_MANY_CONNECTIONS      12	/* 0xc */
#define SSH2_DISCONNECT_AUTH_CANCELLED_BY_USER    13	/* 0xd */
#define SSH2_DISCONNECT_NO_MORE_AUTH_METHODS_AVAILABLE 14	/* 0xe */
#define SSH2_DISCONNECT_ILLEGAL_USER_NAME         15	/* 0xf */

static const char *const ssh2_disconnect_reasons[] = {
    NULL,
    "host not allowed to connect",
    "protocol error",
    "key exchange failed",
    "host authentication failed",
    "MAC error",
    "compression error",
    "service not available",
    "protocol version not supported",
    "host key not verifiable",
    "connection lost",
    "by application",
    "too many connections",
    "auth cancelled by user",
    "no more auth methods available",
    "illegal user name",
};

#define SSH2_OPEN_ADMINISTRATIVELY_PROHIBITED     1	/* 0x1 */
#define SSH2_OPEN_CONNECT_FAILED                  2	/* 0x2 */
#define SSH2_OPEN_UNKNOWN_CHANNEL_TYPE            3	/* 0x3 */
#define SSH2_OPEN_RESOURCE_SHORTAGE               4	/* 0x4 */

#define SSH2_EXTENDED_DATA_STDERR                 1	/* 0x1 */

/*
 * Various remote-bug flags.
 */
#define BUG_CHOKES_ON_SSH1_IGNORE                 1
#define BUG_SSH2_HMAC                             2
#define BUG_NEEDS_SSH1_PLAIN_PASSWORD        	  4
#define BUG_CHOKES_ON_RSA	        	  8
#define BUG_SSH2_RSA_PADDING	        	 16
#define BUG_SSH2_DERIVEKEY                       32
#define BUG_SSH2_REKEY                           64
#define BUG_SSH2_PK_SESSIONID                   128
#define BUG_SSH2_MAXPKT				256
#define BUG_CHOKES_ON_SSH2_IGNORE               512

/*
 * Codes for terminal modes.
 * Most of these are the same in SSH-1 and SSH-2.
 * This list is derived from RFC 4254 and
 * SSH-1 RFC-1.2.31.
 */
static const struct {
    const char* const mode;
    int opcode;
    enum { TTY_OP_CHAR, TTY_OP_BOOL } type;
} ssh_ttymodes[] = {
    /* "V" prefix discarded for special characters relative to SSH specs */
    { "INTR",	      1, TTY_OP_CHAR },
    { "QUIT",	      2, TTY_OP_CHAR },
    { "ERASE",	      3, TTY_OP_CHAR },
    { "KILL",	      4, TTY_OP_CHAR },
    { "EOF",	      5, TTY_OP_CHAR },
    { "EOL",	      6, TTY_OP_CHAR },
    { "EOL2",	      7, TTY_OP_CHAR },
    { "START",	      8, TTY_OP_CHAR },
    { "STOP",	      9, TTY_OP_CHAR },
    { "SUSP",	     10, TTY_OP_CHAR },
    { "DSUSP",	     11, TTY_OP_CHAR },
    { "REPRINT",     12, TTY_OP_CHAR },
    { "WERASE",	     13, TTY_OP_CHAR },
    { "LNEXT",	     14, TTY_OP_CHAR },
    { "FLUSH",	     15, TTY_OP_CHAR },
    { "SWTCH",	     16, TTY_OP_CHAR },
    { "STATUS",	     17, TTY_OP_CHAR },
    { "DISCARD",     18, TTY_OP_CHAR },
    { "IGNPAR",	     30, TTY_OP_BOOL },
    { "PARMRK",	     31, TTY_OP_BOOL },
    { "INPCK",	     32, TTY_OP_BOOL },
    { "ISTRIP",	     33, TTY_OP_BOOL },
    { "INLCR",	     34, TTY_OP_BOOL },
    { "IGNCR",	     35, TTY_OP_BOOL },
    { "ICRNL",	     36, TTY_OP_BOOL },
    { "IUCLC",	     37, TTY_OP_BOOL },
    { "IXON",	     38, TTY_OP_BOOL },
    { "IXANY",	     39, TTY_OP_BOOL },
    { "IXOFF",	     40, TTY_OP_BOOL },
    { "IMAXBEL",     41, TTY_OP_BOOL },
    { "ISIG",	     50, TTY_OP_BOOL },
    { "ICANON",	     51, TTY_OP_BOOL },
    { "XCASE",	     52, TTY_OP_BOOL },
    { "ECHO",	     53, TTY_OP_BOOL },
    { "ECHOE",	     54, TTY_OP_BOOL },
    { "ECHOK",	     55, TTY_OP_BOOL },
    { "ECHONL",	     56, TTY_OP_BOOL },
    { "NOFLSH",	     57, TTY_OP_BOOL },
    { "TOSTOP",	     58, TTY_OP_BOOL },
    { "IEXTEN",	     59, TTY_OP_BOOL },
    { "ECHOCTL",     60, TTY_OP_BOOL },
    { "ECHOKE",	     61, TTY_OP_BOOL },
    { "PENDIN",	     62, TTY_OP_BOOL }, /* XXX is this a real mode? */
    { "OPOST",	     70, TTY_OP_BOOL },
    { "OLCUC",	     71, TTY_OP_BOOL },
    { "ONLCR",	     72, TTY_OP_BOOL },
    { "OCRNL",	     73, TTY_OP_BOOL },
    { "ONOCR",	     74, TTY_OP_BOOL },
    { "ONLRET",	     75, TTY_OP_BOOL },
    { "CS7",	     90, TTY_OP_BOOL },
    { "CS8",	     91, TTY_OP_BOOL },
    { "PARENB",	     92, TTY_OP_BOOL },
    { "PARODD",	     93, TTY_OP_BOOL }
};

/* Miscellaneous other tty-related constants. */
#define SSH_TTY_OP_END		  0
/* The opcodes for ISPEED/OSPEED differ between SSH-1 and SSH-2. */
#define SSH1_TTY_OP_ISPEED	192
#define SSH1_TTY_OP_OSPEED	193
#define SSH2_TTY_OP_ISPEED	128
#define SSH2_TTY_OP_OSPEED	129

/* Helper functions for parsing tty-related config. */
static unsigned int ssh_tty_parse_specchar(char *s)
{
    unsigned int ret;
    if (*s) {
	char *next = NULL;
	ret = ctrlparse(s, &next);
	if (!next) ret = s[0];
    } else {
	ret = 255; /* special value meaning "don't set" */
    }
    return ret;
}
static unsigned int ssh_tty_parse_boolean(char *s)
{
    if (stricmp(s, "yes") == 0 ||
	stricmp(s, "on") == 0 ||
	stricmp(s, "true") == 0 ||
	stricmp(s, "+") == 0)
	return 1; /* true */
    else if (stricmp(s, "no") == 0 ||
	     stricmp(s, "off") == 0 ||
	     stricmp(s, "false") == 0 ||
	     stricmp(s, "-") == 0)
	return 0; /* false */
    else
	return (atoi(s) != 0);
}

#define translate(x) if (type == x) return #x
#define translatek(x,ctx) if (type == x && (pkt_kctx == ctx)) return #x
#define translatea(x,ctx) if (type == x && (pkt_actx == ctx)) return #x
static char *ssh1_pkt_type(int type)
{
    translate(SSH1_MSG_DISCONNECT);
    translate(SSH1_SMSG_PUBLIC_KEY);
    translate(SSH1_CMSG_SESSION_KEY);
    translate(SSH1_CMSG_USER);
    translate(SSH1_CMSG_AUTH_RSA);
    translate(SSH1_SMSG_AUTH_RSA_CHALLENGE);
    translate(SSH1_CMSG_AUTH_RSA_RESPONSE);
    translate(SSH1_CMSG_AUTH_PASSWORD);
    translate(SSH1_CMSG_REQUEST_PTY);
    translate(SSH1_CMSG_WINDOW_SIZE);
    translate(SSH1_CMSG_EXEC_SHELL);
    translate(SSH1_CMSG_EXEC_CMD);
    translate(SSH1_SMSG_SUCCESS);
    translate(SSH1_SMSG_FAILURE);
    translate(SSH1_CMSG_STDIN_DATA);
    translate(SSH1_SMSG_STDOUT_DATA);
    translate(SSH1_SMSG_STDERR_DATA);
    translate(SSH1_CMSG_EOF);
    translate(SSH1_SMSG_EXIT_STATUS);
    translate(SSH1_MSG_CHANNEL_OPEN_CONFIRMATION);
    translate(SSH1_MSG_CHANNEL_OPEN_FAILURE);
    translate(SSH1_MSG_CHANNEL_DATA);
    translate(SSH1_MSG_CHANNEL_CLOSE);
    translate(SSH1_MSG_CHANNEL_CLOSE_CONFIRMATION);
    translate(SSH1_SMSG_X11_OPEN);
    translate(SSH1_CMSG_PORT_FORWARD_REQUEST);
    translate(SSH1_MSG_PORT_OPEN);
    translate(SSH1_CMSG_AGENT_REQUEST_FORWARDING);
    translate(SSH1_SMSG_AGENT_OPEN);
    translate(SSH1_MSG_IGNORE);
    translate(SSH1_CMSG_EXIT_CONFIRMATION);
    translate(SSH1_CMSG_X11_REQUEST_FORWARDING);
    translate(SSH1_CMSG_AUTH_RHOSTS_RSA);
    translate(SSH1_MSG_DEBUG);
    translate(SSH1_CMSG_REQUEST_COMPRESSION);
    translate(SSH1_CMSG_AUTH_TIS);
    translate(SSH1_SMSG_AUTH_TIS_CHALLENGE);
    translate(SSH1_CMSG_AUTH_TIS_RESPONSE);
    translate(SSH1_CMSG_AUTH_CCARD);
    translate(SSH1_SMSG_AUTH_CCARD_CHALLENGE);
    translate(SSH1_CMSG_AUTH_CCARD_RESPONSE);
    return "unknown";
}
static char *ssh2_pkt_type(Pkt_KCtx pkt_kctx, Pkt_ACtx pkt_actx, int type)
{
    translatea(SSH2_MSG_USERAUTH_GSSAPI_RESPONSE,SSH2_PKTCTX_GSSAPI);
    translatea(SSH2_MSG_USERAUTH_GSSAPI_TOKEN,SSH2_PKTCTX_GSSAPI);
    translatea(SSH2_MSG_USERAUTH_GSSAPI_EXCHANGE_COMPLETE,SSH2_PKTCTX_GSSAPI);
    translatea(SSH2_MSG_USERAUTH_GSSAPI_ERROR,SSH2_PKTCTX_GSSAPI);
    translatea(SSH2_MSG_USERAUTH_GSSAPI_ERRTOK,SSH2_PKTCTX_GSSAPI);
    translatea(SSH2_MSG_USERAUTH_GSSAPI_MIC, SSH2_PKTCTX_GSSAPI);
    translate(SSH2_MSG_DISCONNECT);
    translate(SSH2_MSG_IGNORE);
    translate(SSH2_MSG_UNIMPLEMENTED);
    translate(SSH2_MSG_DEBUG);
    translate(SSH2_MSG_SERVICE_REQUEST);
    translate(SSH2_MSG_SERVICE_ACCEPT);
    translate(SSH2_MSG_KEXINIT);
    translate(SSH2_MSG_NEWKEYS);
    translatek(SSH2_MSG_KEXDH_INIT, SSH2_PKTCTX_DHGROUP);
    translatek(SSH2_MSG_KEXDH_REPLY, SSH2_PKTCTX_DHGROUP);
    translatek(SSH2_MSG_KEX_DH_GEX_REQUEST, SSH2_PKTCTX_DHGEX);
    translatek(SSH2_MSG_KEX_DH_GEX_GROUP, SSH2_PKTCTX_DHGEX);
    translatek(SSH2_MSG_KEX_DH_GEX_INIT, SSH2_PKTCTX_DHGEX);
    translatek(SSH2_MSG_KEX_DH_GEX_REPLY, SSH2_PKTCTX_DHGEX);
    translatek(SSH2_MSG_KEXRSA_PUBKEY, SSH2_PKTCTX_RSAKEX);
    translatek(SSH2_MSG_KEXRSA_SECRET, SSH2_PKTCTX_RSAKEX);
    translatek(SSH2_MSG_KEXRSA_DONE, SSH2_PKTCTX_RSAKEX);
    translate(SSH2_MSG_USERAUTH_REQUEST);
    translate(SSH2_MSG_USERAUTH_FAILURE);
    translate(SSH2_MSG_USERAUTH_SUCCESS);
    translate(SSH2_MSG_USERAUTH_BANNER);
    translatea(SSH2_MSG_USERAUTH_PK_OK, SSH2_PKTCTX_PUBLICKEY);
    translatea(SSH2_MSG_USERAUTH_PASSWD_CHANGEREQ, SSH2_PKTCTX_PASSWORD);
    translatea(SSH2_MSG_USERAUTH_INFO_REQUEST, SSH2_PKTCTX_KBDINTER);
    translatea(SSH2_MSG_USERAUTH_INFO_RESPONSE, SSH2_PKTCTX_KBDINTER);
    translate(SSH2_MSG_GLOBAL_REQUEST);
    translate(SSH2_MSG_REQUEST_SUCCESS);
    translate(SSH2_MSG_REQUEST_FAILURE);
    translate(SSH2_MSG_CHANNEL_OPEN);
    translate(SSH2_MSG_CHANNEL_OPEN_CONFIRMATION);
    translate(SSH2_MSG_CHANNEL_OPEN_FAILURE);
    translate(SSH2_MSG_CHANNEL_WINDOW_ADJUST);
    translate(SSH2_MSG_CHANNEL_DATA);
    translate(SSH2_MSG_CHANNEL_EXTENDED_DATA);
    translate(SSH2_MSG_CHANNEL_EOF);
    translate(SSH2_MSG_CHANNEL_CLOSE);
    translate(SSH2_MSG_CHANNEL_REQUEST);
    translate(SSH2_MSG_CHANNEL_SUCCESS);
    translate(SSH2_MSG_CHANNEL_FAILURE);
    return "unknown";
}
#undef translate
#undef translatec

/* Enumeration values for fields in SSH-1 packets */
enum {
    PKT_END, PKT_INT, PKT_CHAR, PKT_DATA, PKT_STR, PKT_BIGNUM,
    /* These values are for communicating relevant semantics of
     * fields to the packet logging code. */
    PKTT_OTHER, PKTT_PASSWORD, PKTT_DATA
};

/*
 * Coroutine mechanics for the sillier bits of the code. If these
 * macros look impenetrable to you, you might find it helpful to
 * read
 * 
 *   http://www.chiark.greenend.org.uk/~sgtatham/coroutines.html
 * 
 * which explains the theory behind these macros.
 * 
 * In particular, if you are getting `case expression not constant'
 * errors when building with MS Visual Studio, this is because MS's
 * Edit and Continue debugging feature causes their compiler to
 * violate ANSI C. To disable Edit and Continue debugging:
 * 
 *  - right-click ssh.c in the FileView
 *  - click Settings
 *  - select the C/C++ tab and the General category
 *  - under `Debug info:', select anything _other_ than `Program
 *    Database for Edit and Continue'.
 */
#define crBegin(v)	{ int *crLine = &v; switch(v) { case 0:;
#define crState(t) \
    struct t *s; \
    if (!ssh->t) ssh->t = snew(struct t); \
    s = ssh->t;
#define crFinish(z)	} *crLine = 0; return (z); }
#define crFinishV	} *crLine = 0; return; }
#define crReturn(z)	\
	do {\
	    *crLine =__LINE__; return (z); case __LINE__:;\
	} while (0)
#define crReturnV	\
	do {\
	    *crLine=__LINE__; return; case __LINE__:;\
	} while (0)
#define crStop(z)	do{ *crLine = 0; return (z); }while(0)
#define crStopV		do{ *crLine = 0; return; }while(0)
#define crWaitUntil(c)	do { crReturn(0); } while (!(c))
#define crWaitUntilV(c)	do { crReturnV; } while (!(c))

typedef struct ssh_tag *Ssh;
struct Packet;

static struct Packet *ssh1_pkt_init(int pkt_type);
static struct Packet *ssh2_pkt_init(int pkt_type);
static void ssh_pkt_ensure(struct Packet *, int length);
static void ssh_pkt_adddata(struct Packet *, void *data, int len);
static void ssh_pkt_addbyte(struct Packet *, unsigned char value);
static void ssh2_pkt_addbool(struct Packet *, unsigned char value);
static void ssh_pkt_adduint32(struct Packet *, unsigned long value);
static void ssh_pkt_addstring_start(struct Packet *);
static void ssh_pkt_addstring_str(struct Packet *, char *data);
static void ssh_pkt_addstring_data(struct Packet *, char *data, int len);
static void ssh_pkt_addstring(struct Packet *, char *data);
static unsigned char *ssh2_mpint_fmt(Bignum b, int *len);
static void ssh1_pkt_addmp(struct Packet *, Bignum b);
static void ssh2_pkt_addmp(struct Packet *, Bignum b);
static int ssh2_pkt_construct(Ssh, struct Packet *);
static void ssh2_pkt_send(Ssh, struct Packet *);
static void ssh2_pkt_send_noqueue(Ssh, struct Packet *);
static int do_ssh1_login(Ssh ssh, unsigned char *in, int inlen,
			 struct Packet *pktin);
static void do_ssh2_authconn(Ssh ssh, unsigned char *in, int inlen,
			     struct Packet *pktin);

/*
 * Buffer management constants. There are several of these for
 * various different purposes:
 * 
 *  - SSH1_BUFFER_LIMIT is the amount of backlog that must build up
 *    on a local data stream before we throttle the whole SSH
 *    connection (in SSH-1 only). Throttling the whole connection is
 *    pretty drastic so we set this high in the hope it won't
 *    happen very often.
 * 
 *  - SSH_MAX_BACKLOG is the amount of backlog that must build up
 *    on the SSH connection itself before we defensively throttle
 *    _all_ local data streams. This is pretty drastic too (though
 *    thankfully unlikely in SSH-2 since the window mechanism should
 *    ensure that the server never has any need to throttle its end
 *    of the connection), so we set this high as well.
 * 
 *  - OUR_V2_WINSIZE is the maximum window size we present on SSH-2
 *    channels.
 *
 *  - OUR_V2_BIGWIN is the window size we advertise for the only
 *    channel in a simple connection.  It must be <= INT_MAX.
 *
 *  - OUR_V2_MAXPKT is the official "maximum packet size" we send
 *    to the remote side. This actually has nothing to do with the
 *    size of the _packet_, but is instead a limit on the amount
 *    of data we're willing to receive in a single SSH2 channel
 *    data message.
 *
 *  - OUR_V2_PACKETLIMIT is actually the maximum size of SSH
 *    _packet_ we're prepared to cope with.  It must be a multiple
 *    of the cipher block size, and must be at least 35000.
 */

#define SSH1_BUFFER_LIMIT 32768
#define SSH_MAX_BACKLOG 32768
#define OUR_V2_WINSIZE 16384
#define OUR_V2_BIGWIN 0x7fffffff
#define OUR_V2_MAXPKT 0x4000UL
#define OUR_V2_PACKETLIMIT 0x9000UL

/* Maximum length of passwords/passphrases (arbitrary) */
#define SSH_MAX_PASSWORD_LEN 100

const static struct ssh_signkey *hostkey_algs[] = { &ssh_rsa, &ssh_dss };

const static struct ssh_mac *macs[] = {
    &ssh_hmac_sha1, &ssh_hmac_sha1_96, &ssh_hmac_md5
};
const static struct ssh_mac *buggymacs[] = {
    &ssh_hmac_sha1_buggy, &ssh_hmac_sha1_96_buggy, &ssh_hmac_md5
};

static void *ssh_comp_none_init(void)
{
    return NULL;
}
static void ssh_comp_none_cleanup(void *handle)
{
}
static int ssh_comp_none_block(void *handle, unsigned char *block, int len,
			       unsigned char **outblock, int *outlen)
{
    return 0;
}
static int ssh_comp_none_disable(void *handle)
{
    return 0;
}
const static struct ssh_compress ssh_comp_none = {
    "none", NULL,
    ssh_comp_none_init, ssh_comp_none_cleanup, ssh_comp_none_block,
    ssh_comp_none_init, ssh_comp_none_cleanup, ssh_comp_none_block,
    ssh_comp_none_disable, NULL
};
extern const struct ssh_compress ssh_zlib;
const static struct ssh_compress *compressions[] = {
    &ssh_zlib, &ssh_comp_none
};

enum {				       /* channel types */
    CHAN_MAINSESSION,
    CHAN_X11,
    CHAN_AGENT,
    CHAN_SOCKDATA,
    CHAN_SOCKDATA_DORMANT	       /* one the remote hasn't confirmed */
};

/*
 * little structure to keep track of outstanding WINDOW_ADJUSTs
 */
struct winadj {
    struct winadj *next;
    unsigned size;
};

/*
 * 2-3-4 tree storing channels.
 */
struct ssh_channel {
    Ssh ssh;			       /* pointer back to main context */
    unsigned remoteid, localid;
    int type;
    /* True if we opened this channel but server hasn't confirmed. */
    int halfopen;
    /*
     * In SSH-1, this value contains four bits:
     * 
     *   1   We have sent SSH1_MSG_CHANNEL_CLOSE.
     *   2   We have sent SSH1_MSG_CHANNEL_CLOSE_CONFIRMATION.
     *   4   We have received SSH1_MSG_CHANNEL_CLOSE.
     *   8   We have received SSH1_MSG_CHANNEL_CLOSE_CONFIRMATION.
     * 
     * A channel is completely finished with when all four bits are set.
     */
    int closes;

    /*
     * This flag indicates that a close is pending on the outgoing
     * side of the channel: that is, wherever we're getting the data
     * for this channel has sent us some data followed by EOF. We
     * can't actually close the channel until we've finished sending
     * the data, so we set this flag instead to remind us to
     * initiate the closing process once our buffer is clear.
     */
    int pending_close;

    /*
     * True if this channel is causing the underlying connection to be
     * throttled.
     */
    int throttling_conn;
    union {
	struct ssh2_data_channel {
	    bufchain outbuffer;
	    unsigned remwindow, remmaxpkt;
	    /* locwindow is signed so we can cope with excess data. */
	    int locwindow, locmaxwin;
	    /*
	     * remlocwin is the amount of local window that we think
	     * the remote end had available to it after it sent the
	     * last data packet or window adjust ack.
	     */
	    int remlocwin;
	    /*
	     * These store the list of window adjusts that haven't
	     * been acked.
	     */
	    struct winadj *winadj_head, *winadj_tail;
	    enum { THROTTLED, UNTHROTTLING, UNTHROTTLED } throttle_state;
	} v2;
    } v;
    union {
	struct ssh_agent_channel {
	    unsigned char *message;
	    unsigned char msglen[4];
	    unsigned lensofar, totallen;
	} a;
	struct ssh_x11_channel {
	    Socket s;
	} x11;
	struct ssh_pfd_channel {
	    Socket s;
	} pfd;
    } u;
};

/*
 * 2-3-4 tree storing remote->local port forwardings. SSH-1 and SSH-2
 * use this structure in different ways, reflecting SSH-2's
 * altogether saner approach to port forwarding.
 * 
 * In SSH-1, you arrange a remote forwarding by sending the server
 * the remote port number, and the local destination host:port.
 * When a connection comes in, the server sends you back that
 * host:port pair, and you connect to it. This is a ready-made
 * security hole if you're not on the ball: a malicious server
 * could send you back _any_ host:port pair, so if you trustingly
 * connect to the address it gives you then you've just opened the
 * entire inside of your corporate network just by connecting
 * through it to a dodgy SSH server. Hence, we must store a list of
 * host:port pairs we _are_ trying to forward to, and reject a
 * connection request from the server if it's not in the list.
 * 
 * In SSH-2, each side of the connection minds its own business and
 * doesn't send unnecessary information to the other. You arrange a
 * remote forwarding by sending the server just the remote port
 * number. When a connection comes in, the server tells you which
 * of its ports was connected to; and _you_ have to remember what
 * local host:port pair went with that port number.
 * 
 * Hence, in SSH-1 this structure is indexed by destination
 * host:port pair, whereas in SSH-2 it is indexed by source port.
 */
struct ssh_portfwd; /* forward declaration */

struct ssh_rportfwd {
    unsigned sport, dport;
    char dhost[256];
    char *sportdesc;
    struct ssh_portfwd *pfrec;
};
#define free_rportfwd(pf) ( \
    ((pf) ? (sfree((pf)->sportdesc)) : (void)0 ), sfree(pf) )

/*
 * Separately to the rportfwd tree (which is for looking up port
 * open requests from the server), a tree of _these_ structures is
 * used to keep track of all the currently open port forwardings,
 * so that we can reconfigure in mid-session if the user requests
 * it.
 */
struct ssh_portfwd {
    enum { DESTROY, KEEP, CREATE } status;
    int type;
    unsigned sport, dport;
    char *saddr, *daddr;
    char *sserv, *dserv;
    struct ssh_rportfwd *remote;
    int addressfamily;
    void *local;
};
#define free_portfwd(pf) ( \
    ((pf) ? (sfree((pf)->saddr), sfree((pf)->daddr), \
	     sfree((pf)->sserv), sfree((pf)->dserv)) : (void)0 ), sfree(pf) )

struct Packet {
    long length;	    /* length of `data' actually used */
    long forcepad;	    /* SSH-2: force padding to at least this length */
    int type;		    /* only used for incoming packets */
    unsigned long sequence; /* SSH-2 incoming sequence number */
    unsigned char *data;    /* allocated storage */
    unsigned char *body;    /* offset of payload within `data' */
    long savedpos;	    /* temporary index into `data' (for strings) */
    long maxlen;	    /* amount of storage allocated for `data' */
    long encrypted_len;	    /* for SSH-2 total-size counting */

    /*
     * State associated with packet logging
     */
    int logmode;
    int nblanks;
    struct logblank_t *blanks;
};

static void ssh1_protocol(Ssh ssh, void *vin, int inlen,
			  struct Packet *pktin);
static void ssh2_protocol(Ssh ssh, void *vin, int inlen,
			  struct Packet *pktin);
static void ssh1_protocol_setup(Ssh ssh);
static void ssh2_protocol_setup(Ssh ssh);
static void ssh_size(void *handle, int width, int height);
static void ssh_special(void *handle, Telnet_Special);
static int ssh2_try_send(struct ssh_channel *c);
static void ssh2_add_channel_data(struct ssh_channel *c, char *buf, int len);
static void ssh_throttle_all(Ssh ssh, int enable, int bufsize);
static void ssh2_set_window(struct ssh_channel *c, int newwin);
static int ssh_sendbuffer(void *handle);
static int ssh_do_close(Ssh ssh, int notify_exit);
static unsigned long ssh_pkt_getuint32(struct Packet *pkt);
static int ssh2_pkt_getbool(struct Packet *pkt);
static void ssh_pkt_getstring(struct Packet *pkt, char **p, int *length);
static void ssh2_timer(void *ctx, long now);
static int do_ssh2_transport(Ssh ssh, void *vin, int inlen,
			     struct Packet *pktin);

struct rdpkt1_state_tag {
    long len, pad, biglen, to_read;
    unsigned long realcrc, gotcrc;
    unsigned char *p;
    int i;
    int chunk;
    struct Packet *pktin;
};

struct rdpkt2_state_tag {
    long len, pad, payload, packetlen, maclen;
    int i;
    int cipherblk;
    unsigned long incoming_sequence;
    struct Packet *pktin;
};

typedef void (*handler_fn_t)(Ssh ssh, struct Packet *pktin);
typedef void (*chandler_fn_t)(Ssh ssh, struct Packet *pktin, void *ctx);

struct queued_handler;
struct queued_handler {
    int msg1, msg2;
    chandler_fn_t handler;
    void *ctx;
    struct queued_handler *next;
};

struct ssh_tag {
    const struct plug_function_table *fn;
    /* the above field _must_ be first in the structure */

    char *v_c, *v_s;
    void *exhash;

    Socket s;

    void *ldisc;
    void *logctx;

    unsigned char session_key[32];
    int v1_compressing;
    int v1_remote_protoflags;
    int v1_local_protoflags;
    int agentfwd_enabled;
    int X11_fwd_enabled;
    int remote_bugs;
    const struct ssh_cipher *cipher;
    void *v1_cipher_ctx;
    void *crcda_ctx;
    const struct ssh2_cipher *cscipher, *sccipher;
    void *cs_cipher_ctx, *sc_cipher_ctx;
    const struct ssh_mac *csmac, *scmac;
    void *cs_mac_ctx, *sc_mac_ctx;
    const struct ssh_compress *cscomp, *sccomp;
    void *cs_comp_ctx, *sc_comp_ctx;
    const struct ssh_kex *kex;
    const struct ssh_signkey *hostkey;
    unsigned char v2_session_id[SSH2_KEX_MAX_HASH_LEN];
    int v2_session_id_len;
    void *kex_ctx;

    char *savedhost;
    int savedport;
    int send_ok;
    int echoing, editing;

    void *frontend;

    int ospeed, ispeed;		       /* temporaries */
    int term_width, term_height;

    tree234 *channels;		       /* indexed by local id */
    struct ssh_channel *mainchan;      /* primary session channel */
    int ncmode;			       /* is primary channel direct-tcpip? */
    int exitcode;
    int close_expected;
    int clean_exit;

    tree234 *rportfwds, *portfwds;

    enum {
	SSH_STATE_PREPACKET,
	SSH_STATE_BEFORE_SIZE,
	SSH_STATE_INTERMED,
	SSH_STATE_SESSION,
	SSH_STATE_CLOSED
    } state;

    int size_needed, eof_needed;

    struct Packet **queue;
    int queuelen, queuesize;
    int queueing;
    unsigned char *deferred_send_data;
    int deferred_len, deferred_size;

    /*
     * Gross hack: pscp will try to start SFTP but fall back to
     * scp1 if that fails. This variable is the means by which
     * scp.c can reach into the SSH code and find out which one it
     * got.
     */
    int fallback_cmd;

    bufchain banner;	/* accumulates banners during do_ssh2_authconn */

    Pkt_KCtx pkt_kctx;
    Pkt_ACtx pkt_actx;

    struct X11Display *x11disp;

    int version;
    int conn_throttle_count;
    int overall_bufsize;
    int throttled_all;
    int v1_stdout_throttling;
    unsigned long v2_outgoing_sequence;

    int ssh1_rdpkt_crstate;
    int ssh2_rdpkt_crstate;
    int do_ssh_init_crstate;
    int ssh_gotdata_crstate;
    int do_ssh1_login_crstate;
    int do_ssh1_connection_crstate;
    int do_ssh2_transport_crstate;
    int do_ssh2_authconn_crstate;

    void *do_ssh_init_state;
    void *do_ssh1_login_state;
    void *do_ssh2_transport_state;
    void *do_ssh2_authconn_state;

    struct rdpkt1_state_tag rdpkt1_state;
    struct rdpkt2_state_tag rdpkt2_state;

    /* SSH-1 and SSH-2 use this for different things, but both use it */
    int protocol_initial_phase_done;

    void (*protocol) (Ssh ssh, void *vin, int inlen,
		      struct Packet *pkt);
    struct Packet *(*s_rdpkt) (Ssh ssh, unsigned char **data, int *datalen);

    /*
     * We maintain a full _copy_ of a Config structure here, not
     * merely a pointer to it. That way, when we're passed a new
     * one for reconfiguration, we can check the differences and
     * potentially reconfigure port forwardings etc in mid-session.
     */
    Config cfg;

    /*
     * Used to transfer data back from async callbacks.
     */
    void *agent_response;
    int agent_response_len;
    int user_response;

    /*
     * The SSH connection can be set as `frozen', meaning we are
     * not currently accepting incoming data from the network. This
     * is slightly more serious than setting the _socket_ as
     * frozen, because we may already have had data passed to us
     * from the network which we need to delay processing until
     * after the freeze is lifted, so we also need a bufchain to
     * store that data.
     */
    int frozen;
    bufchain queued_incoming_data;

    /*
     * Dispatch table for packet types that we may have to deal
     * with at any time.
     */
    handler_fn_t packet_dispatch[256];

    /*
     * Queues of one-off handler functions for success/failure
     * indications from a request.
     */
    struct queued_handler *qhead, *qtail;

    /*
     * This module deals with sending keepalives.
     */
    Pinger pinger;

    /*
     * Track incoming and outgoing data sizes and time, for
     * size-based rekeys.
     */
    unsigned long incoming_data_size, outgoing_data_size, deferred_data_size;
    unsigned long max_data_size;
    int kex_in_progress;
    long next_rekey, last_rekey;
    char *deferred_rekey_reason;    /* points to STATIC string; don't free */

    /*
     * Fully qualified host name, which we need if doing GSSAPI.
     */
    char *fullhostname;

#ifndef NO_GSSAPI
    /*
     * GSSAPI libraries for this session.
     */
    struct ssh_gss_liblist *gsslibs;
#endif
};

#define logevent(s) logevent(ssh->frontend, s)

/* logevent, only printf-formatted. */
static void logeventf(Ssh ssh, const char *fmt, ...)
{
    va_list ap;
    char *buf;

    va_start(ap, fmt);
    buf = dupvprintf(fmt, ap);
    va_end(ap);
    logevent(buf);
    sfree(buf);
}

#define bombout(msg) \
    do { \
        char *text = dupprintf msg; \
	ssh_do_close(ssh, FALSE); \
        logevent(text); \
        connection_fatal(ssh->frontend, "%s", text); \
        sfree(text); \
    } while (0)

/* Functions to leave bits out of the SSH packet log file. */

static void dont_log_password(Ssh ssh, struct Packet *pkt, int blanktype)
{
    if (ssh->cfg.logomitpass)
	pkt->logmode = blanktype;
}

static void dont_log_data(Ssh ssh, struct Packet *pkt, int blanktype)
{
    if (ssh->cfg.logomitdata)
	pkt->logmode = blanktype;
}

static void end_log_omission(Ssh ssh, struct Packet *pkt)
{
    pkt->logmode = PKTLOG_EMIT;
}

/* Helper function for common bits of parsing cfg.ttymodes. */
static void parse_ttymodes(Ssh ssh, char *modes,
			   void (*do_mode)(void *data, char *mode, char *val),
			   void *data)
{
    while (*modes) {
	char *t = strchr(modes, '\t');
	char *m = snewn(t-modes+1, char);
	char *val;
	strncpy(m, modes, t-modes);
	m[t-modes] = '\0';
	if (*(t+1) == 'A')
	    val = get_ttymode(ssh->frontend, m);
	else
	    val = dupstr(t+2);
	if (val)
	    do_mode(data, m, val);
	sfree(m);
	sfree(val);
	modes += strlen(modes) + 1;
    }
}

static int ssh_channelcmp(void *av, void *bv)
{
    struct ssh_channel *a = (struct ssh_channel *) av;
    struct ssh_channel *b = (struct ssh_channel *) bv;
    if (a->localid < b->localid)
	return -1;
    if (a->localid > b->localid)
	return +1;
    return 0;
}
static int ssh_channelfind(void *av, void *bv)
{
    unsigned *a = (unsigned *) av;
    struct ssh_channel *b = (struct ssh_channel *) bv;
    if (*a < b->localid)
	return -1;
    if (*a > b->localid)
	return +1;
    return 0;
}

static int ssh_rportcmp_ssh1(void *av, void *bv)
{
    struct ssh_rportfwd *a = (struct ssh_rportfwd *) av;
    struct ssh_rportfwd *b = (struct ssh_rportfwd *) bv;
    int i;
    if ( (i = strcmp(a->dhost, b->dhost)) != 0)
	return i < 0 ? -1 : +1;
    if (a->dport > b->dport)
	return +1;
    if (a->dport < b->dport)
	return -1;
    return 0;
}

static int ssh_rportcmp_ssh2(void *av, void *bv)
{
    struct ssh_rportfwd *a = (struct ssh_rportfwd *) av;
    struct ssh_rportfwd *b = (struct ssh_rportfwd *) bv;

    if (a->sport > b->sport)
	return +1;
    if (a->sport < b->sport)
	return -1;
    return 0;
}

/*
 * Special form of strcmp which can cope with NULL inputs. NULL is
 * defined to sort before even the empty string.
 */
static int nullstrcmp(const char *a, const char *b)
{
    if (a == NULL && b == NULL)
	return 0;
    if (a == NULL)
	return -1;
    if (b == NULL)
	return +1;
    return strcmp(a, b);
}

static int ssh_portcmp(void *av, void *bv)
{
    struct ssh_portfwd *a = (struct ssh_portfwd *) av;
    struct ssh_portfwd *b = (struct ssh_portfwd *) bv;
    int i;
    if (a->type > b->type)
	return +1;
    if (a->type < b->type)
	return -1;
    if (a->addressfamily > b->addressfamily)
	return +1;
    if (a->addressfamily < b->addressfamily)
	return -1;
    if ( (i = nullstrcmp(a->saddr, b->saddr)) != 0)
	return i < 0 ? -1 : +1;
    if (a->sport > b->sport)
	return +1;
    if (a->sport < b->sport)
	return -1;
    if (a->type != 'D') {
	if ( (i = nullstrcmp(a->daddr, b->daddr)) != 0)
	    return i < 0 ? -1 : +1;
	if (a->dport > b->dport)
	    return +1;
	if (a->dport < b->dport)
	    return -1;
    }
    return 0;
}

static int alloc_channel_id(Ssh ssh)
{
    const unsigned CHANNEL_NUMBER_OFFSET = 256;
    unsigned low, high, mid;
    int tsize;
    struct ssh_channel *c;

    /*
     * First-fit allocation of channel numbers: always pick the
     * lowest unused one. To do this, binary-search using the
     * counted B-tree to find the largest channel ID which is in a
     * contiguous sequence from the beginning. (Precisely
     * everything in that sequence must have ID equal to its tree
     * index plus CHANNEL_NUMBER_OFFSET.)
     */
    tsize = count234(ssh->channels);

    low = -1;
    high = tsize;
    while (high - low > 1) {
	mid = (high + low) / 2;
	c = index234(ssh->channels, mid);
	if (c->localid == mid + CHANNEL_NUMBER_OFFSET)
	    low = mid;		       /* this one is fine */
	else
	    high = mid;		       /* this one is past it */
    }
    /*
     * Now low points to either -1, or the tree index of the
     * largest ID in the initial sequence.
     */
    {
	unsigned i = low + 1 + CHANNEL_NUMBER_OFFSET;
	assert(NULL == find234(ssh->channels, &i, ssh_channelfind));
    }
    return low + 1 + CHANNEL_NUMBER_OFFSET;
}

static void c_write_stderr(int trusted, const char *buf, int len)
{
    int i;
    for (i = 0; i < len; i++)
	if (buf[i] != '\r' && (trusted || buf[i] == '\n' || (buf[i] & 0x60)))
	    fputc(buf[i], stderr);
}

static void c_write(Ssh ssh, const char *buf, int len)
{
    if (flags & FLAG_STDERR)
	c_write_stderr(1, buf, len);
    else
	from_backend(ssh->frontend, 1, buf, len);
}

static void c_write_untrusted(Ssh ssh, const char *buf, int len)
{
    if (flags & FLAG_STDERR)
	c_write_stderr(0, buf, len);
    else
	from_backend_untrusted(ssh->frontend, buf, len);
}

static void c_write_str(Ssh ssh, const char *buf)
{
    c_write(ssh, buf, strlen(buf));
}

static void ssh_free_packet(struct Packet *pkt)
{
    sfree(pkt->data);
    sfree(pkt);
}
static struct Packet *ssh_new_packet(void)
{
    struct Packet *pkt = snew(struct Packet);

    pkt->body = pkt->data = NULL;
    pkt->maxlen = 0;
    pkt->logmode = PKTLOG_EMIT;
    pkt->nblanks = 0;
    pkt->blanks = NULL;

    return pkt;
}

/*
 * Collect incoming data in the incoming packet buffer.
 * Decipher and verify the packet when it is completely read.
 * Drop SSH1_MSG_DEBUG and SSH1_MSG_IGNORE packets.
 * Update the *data and *datalen variables.
 * Return a Packet structure when a packet is completed.
 */
static struct Packet *ssh1_rdpkt(Ssh ssh, unsigned char **data, int *datalen)
{
    struct rdpkt1_state_tag *st = &ssh->rdpkt1_state;

    crBegin(ssh->ssh1_rdpkt_crstate);

    st->pktin = ssh_new_packet();

    st->pktin->type = 0;
    st->pktin->length = 0;

    for (st->i = st->len = 0; st->i < 4; st->i++) {
	while ((*datalen) == 0)
	    crReturn(NULL);
	st->len = (st->len << 8) + **data;
	(*data)++, (*datalen)--;
    }

    st->pad = 8 - (st->len % 8);
    st->biglen = st->len + st->pad;
    st->pktin->length = st->len - 5;

    if (st->biglen < 0) {
        bombout(("Extremely large packet length from server suggests"
		 " data stream corruption"));
	ssh_free_packet(st->pktin);
        crStop(NULL);
    }

    st->pktin->maxlen = st->biglen;
    st->pktin->data = snewn(st->biglen + APIEXTRA, unsigned char);

    st->to_read = st->biglen;
    st->p = st->pktin->data;
    while (st->to_read > 0) {
	st->chunk = st->to_read;
	while ((*datalen) == 0)
	    crReturn(NULL);
	if (st->chunk > (*datalen))
	    st->chunk = (*datalen);
	memcpy(st->p, *data, st->chunk);
	*data += st->chunk;
	*datalen -= st->chunk;
	st->p += st->chunk;
	st->to_read -= st->chunk;
    }

    if (ssh->cipher && detect_attack(ssh->crcda_ctx, st->pktin->data,
				     st->biglen, NULL)) {
        bombout(("Network attack (CRC compensation) detected!"));
	ssh_free_packet(st->pktin);
        crStop(NULL);
    }

    if (ssh->cipher)
	ssh->cipher->decrypt(ssh->v1_cipher_ctx, st->pktin->data, st->biglen);

    st->realcrc = crc32_compute(st->pktin->data, st->biglen - 4);
    st->gotcrc = GET_32BIT(st->pktin->data + st->biglen - 4);
    if (st->gotcrc != st->realcrc) {
	bombout(("Incorrect CRC received on packet"));
	ssh_free_packet(st->pktin);
	crStop(NULL);
    }

    st->pktin->body = st->pktin->data + st->pad + 1;
    st->pktin->savedpos = 0;

    if (ssh->v1_compressing) {
	unsigned char *decompblk;
	int decomplen;
	if (!zlib_decompress_block(ssh->sc_comp_ctx,
				   st->pktin->body - 1, st->pktin->length + 1,
				   &decompblk, &decomplen)) {
	    bombout(("Zlib decompression encountered invalid data"));
	    ssh_free_packet(st->pktin);
	    crStop(NULL);
	}

	if (st->pktin->maxlen < st->pad + decomplen) {
	    st->pktin->maxlen = st->pad + decomplen;
	    st->pktin->data = sresize(st->pktin->data,
				      st->pktin->maxlen + APIEXTRA,
				      unsigned char);
	    st->pktin->body = st->pktin->data + st->pad + 1;
	}

	memcpy(st->pktin->body - 1, decompblk, decomplen);
	sfree(decompblk);
	st->pktin->length = decomplen - 1;
    }

    st->pktin->type = st->pktin->body[-1];

    /*
     * Log incoming packet, possibly omitting sensitive fields.
     */
    if (ssh->logctx) {
	int nblanks = 0;
	struct logblank_t blank;
	if (ssh->cfg.logomitdata) {
	    int do_blank = FALSE, blank_prefix = 0;
	    /* "Session data" packets - omit the data field */
	    if ((st->pktin->type == SSH1_SMSG_STDOUT_DATA) ||
		(st->pktin->type == SSH1_SMSG_STDERR_DATA)) {
		do_blank = TRUE; blank_prefix = 4;
	    } else if (st->pktin->type == SSH1_MSG_CHANNEL_DATA) {
		do_blank = TRUE; blank_prefix = 8;
	    }
	    if (do_blank) {
		blank.offset = blank_prefix;
		blank.len = st->pktin->length;
		blank.type = PKTLOG_OMIT;
		nblanks = 1;
	    }
	}
	log_packet(ssh->logctx,
		   PKT_INCOMING, st->pktin->type,
		   ssh1_pkt_type(st->pktin->type),
		   st->pktin->body, st->pktin->length,
		   nblanks, &blank, NULL);
    }

    crFinish(st->pktin);
}

static struct Packet *ssh2_rdpkt(Ssh ssh, unsigned char **data, int *datalen)
{
    struct rdpkt2_state_tag *st = &ssh->rdpkt2_state;

    crBegin(ssh->ssh2_rdpkt_crstate);

    st->pktin = ssh_new_packet();

    st->pktin->type = 0;
    st->pktin->length = 0;
    if (ssh->sccipher)
	st->cipherblk = ssh->sccipher->blksize;
    else
	st->cipherblk = 8;
    if (st->cipherblk < 8)
	st->cipherblk = 8;
    st->maclen = ssh->scmac ? ssh->scmac->len : 0;

    if (ssh->sccipher && (ssh->sccipher->flags & SSH_CIPHER_IS_CBC) &&
	ssh->scmac) {
	/*
	 * When dealing with a CBC-mode cipher, we want to avoid the
	 * possibility of an attacker's tweaking the ciphertext stream
	 * so as to cause us to feed the same block to the block
	 * cipher more than once and thus leak information
	 * (VU#958563).  The way we do this is not to take any
	 * decisions on the basis of anything we've decrypted until
	 * we've verified it with a MAC.  That includes the packet
	 * length, so we just read data and check the MAC repeatedly,
	 * and when the MAC passes, see if the length we've got is
	 * plausible.
	 */

	/* May as well allocate the whole lot now. */
	st->pktin->data = snewn(OUR_V2_PACKETLIMIT + st->maclen + APIEXTRA,
				unsigned char);

	/* Read an amount corresponding to the MAC. */
	for (st->i = 0; st->i < st->maclen; st->i++) {
	    while ((*datalen) == 0)
		crReturn(NULL);
	    st->pktin->data[st->i] = *(*data)++;
	    (*datalen)--;
	}

	st->packetlen = 0;
	{
	    unsigned char seq[4];
	    ssh->scmac->start(ssh->sc_mac_ctx);
	    PUT_32BIT(seq, st->incoming_sequence);
	    ssh->scmac->bytes(ssh->sc_mac_ctx, seq, 4);
	}

	for (;;) { /* Once around this loop per cipher block. */
	    /* Read another cipher-block's worth, and tack it onto the end. */
	    for (st->i = 0; st->i < st->cipherblk; st->i++) {
		while ((*datalen) == 0)
		    crReturn(NULL);
		st->pktin->data[st->packetlen+st->maclen+st->i] = *(*data)++;
		(*datalen)--;
	    }
	    /* Decrypt one more block (a little further back in the stream). */
	    ssh->sccipher->decrypt(ssh->sc_cipher_ctx,
				   st->pktin->data + st->packetlen,
				   st->cipherblk);
	    /* Feed that block to the MAC. */
	    ssh->scmac->bytes(ssh->sc_mac_ctx,
			      st->pktin->data + st->packetlen, st->cipherblk);
	    st->packetlen += st->cipherblk;
	    /* See if that gives us a valid packet. */
	    if (ssh->scmac->verresult(ssh->sc_mac_ctx,
				      st->pktin->data + st->packetlen) &&
		(st->len = GET_32BIT(st->pktin->data)) + 4 == st->packetlen)
		    break;
	    if (st->packetlen >= OUR_V2_PACKETLIMIT) {
		bombout(("No valid incoming packet found"));
		ssh_free_packet(st->pktin);
		crStop(NULL);
	    }	    
	}
	st->pktin->maxlen = st->packetlen + st->maclen;
	st->pktin->data = sresize(st->pktin->data,
				  st->pktin->maxlen + APIEXTRA,
				  unsigned char);
    } else {
	st->pktin->data = snewn(st->cipherblk + APIEXTRA, unsigned char);

	/*
	 * Acquire and decrypt the first block of the packet. This will
	 * contain the length and padding details.
	 */
	for (st->i = st->len = 0; st->i < st->cipherblk; st->i++) {
	    while ((*datalen) == 0)
		crReturn(NULL);
	    st->pktin->data[st->i] = *(*data)++;
	    (*datalen)--;
	}

	if (ssh->sccipher)
	    ssh->sccipher->decrypt(ssh->sc_cipher_ctx,
				   st->pktin->data, st->cipherblk);

	/*
	 * Now get the length figure.
	 */
	st->len = GET_32BIT(st->pktin->data);

	/*
	 * _Completely_ silly lengths should be stomped on before they
	 * do us any more damage.
	 */
	if (st->len < 0 || st->len > OUR_V2_PACKETLIMIT ||
	    (st->len + 4) % st->cipherblk != 0) {
	    bombout(("Incoming packet was garbled on decryption"));
	    ssh_free_packet(st->pktin);
	    crStop(NULL);
	}

	/*
	 * So now we can work out the total packet length.
	 */
	st->packetlen = st->len + 4;

	/*
	 * Allocate memory for the rest of the packet.
	 */
	st->pktin->maxlen = st->packetlen + st->maclen;
	st->pktin->data = sresize(st->pktin->data,
				  st->pktin->maxlen + APIEXTRA,
				  unsigned char);

	/*
	 * Read and decrypt the remainder of the packet.
	 */
	for (st->i = st->cipherblk; st->i < st->packetlen + st->maclen;
	     st->i++) {
	    while ((*datalen) == 0)
		crReturn(NULL);
	    st->pktin->data[st->i] = *(*data)++;
	    (*datalen)--;
	}
	/* Decrypt everything _except_ the MAC. */
	if (ssh->sccipher)
	    ssh->sccipher->decrypt(ssh->sc_cipher_ctx,
				   st->pktin->data + st->cipherblk,
				   st->packetlen - st->cipherblk);

	/*
	 * Check the MAC.
	 */
	if (ssh->scmac
	    && !ssh->scmac->verify(ssh->sc_mac_ctx, st->pktin->data,
				   st->len + 4, st->incoming_sequence)) {
	    bombout(("Incorrect MAC received on packet"));
	    ssh_free_packet(st->pktin);
	    crStop(NULL);
	}
    }
    /* Get and sanity-check the amount of random padding. */
    st->pad = st->pktin->data[4];
    if (st->pad < 4 || st->len - st->pad < 1) {
	bombout(("Invalid padding length on received packet"));
	ssh_free_packet(st->pktin);
	crStop(NULL);
    }
    /*
     * This enables us to deduce the payload length.
     */
    st->payload = st->len - st->pad - 1;

    st->pktin->length = st->payload + 5;
    st->pktin->encrypted_len = st->packetlen;

    st->pktin->sequence = st->incoming_sequence++;

    /*
     * Decompress packet payload.
     */
    {
	unsigned char *newpayload;
	int newlen;
	if (ssh->sccomp &&
	    ssh->sccomp->decompress(ssh->sc_comp_ctx,
				    st->pktin->data + 5, st->pktin->length - 5,
				    &newpayload, &newlen)) {
	    if (st->pktin->maxlen < newlen + 5) {
		st->pktin->maxlen = newlen + 5;
		st->pktin->data = sresize(st->pktin->data,
					  st->pktin->maxlen + APIEXTRA,
					  unsigned char);
	    }
	    st->pktin->length = 5 + newlen;
	    memcpy(st->pktin->data + 5, newpayload, newlen);
	    sfree(newpayload);
	}
    }

    st->pktin->savedpos = 6;
    st->pktin->body = st->pktin->data;
    st->pktin->type = st->pktin->data[5];

    /*
     * Log incoming packet, possibly omitting sensitive fields.
     */
    if (ssh->logctx) {
	int nblanks = 0;
	struct logblank_t blank;
	if (ssh->cfg.logomitdata) {
	    int do_blank = FALSE, blank_prefix = 0;
	    /* "Session data" packets - omit the data field */
	    if (st->pktin->type == SSH2_MSG_CHANNEL_DATA) {
		do_blank = TRUE; blank_prefix = 8;
	    } else if (st->pktin->type == SSH2_MSG_CHANNEL_EXTENDED_DATA) {
		do_blank = TRUE; blank_prefix = 12;
	    }
	    if (do_blank) {
		blank.offset = blank_prefix;
		blank.len = (st->pktin->length-6) - blank_prefix;
		blank.type = PKTLOG_OMIT;
		nblanks = 1;
	    }
	}
	log_packet(ssh->logctx, PKT_INCOMING, st->pktin->type,
		   ssh2_pkt_type(ssh->pkt_kctx, ssh->pkt_actx,
				 st->pktin->type),
		   st->pktin->data+6, st->pktin->length-6,
		   nblanks, &blank, &st->pktin->sequence);
    }

    crFinish(st->pktin);
}

static int s_wrpkt_prepare(Ssh ssh, struct Packet *pkt, int *offset_p)
{
    int pad, biglen, i, pktoffs;
    unsigned long crc;
#ifdef __SC__
    /*
     * XXX various versions of SC (including 8.8.4) screw up the
     * register allocation in this function and use the same register
     * (D6) for len and as a temporary, with predictable results.  The
     * following sledgehammer prevents this.
     */
    volatile
#endif
    int len;

    if (ssh->logctx)
	log_packet(ssh->logctx, PKT_OUTGOING, pkt->data[12],
		   ssh1_pkt_type(pkt->data[12]),
		   pkt->body, pkt->length - (pkt->body - pkt->data),
		   pkt->nblanks, pkt->blanks, NULL);
    sfree(pkt->blanks); pkt->blanks = NULL;
    pkt->nblanks = 0;

    if (ssh->v1_compressing) {
	unsigned char *compblk;
	int complen;
	zlib_compress_block(ssh->cs_comp_ctx,
			    pkt->data + 12, pkt->length - 12,
			    &compblk, &complen);
	ssh_pkt_ensure(pkt, complen + 2);   /* just in case it's got bigger */
	memcpy(pkt->data + 12, compblk, complen);
	sfree(compblk);
	pkt->length = complen + 12;
    }

    ssh_pkt_ensure(pkt, pkt->length + 4); /* space for CRC */
    pkt->length += 4;
    len = pkt->length - 4 - 8;	/* len(type+data+CRC) */
    pad = 8 - (len % 8);
    pktoffs = 8 - pad;
    biglen = len + pad;		/* len(padding+type+data+CRC) */

    for (i = pktoffs; i < 4+8; i++)
	pkt->data[i] = random_byte();
    crc = crc32_compute(pkt->data + pktoffs + 4, biglen - 4); /* all ex len */
    PUT_32BIT(pkt->data + pktoffs + 4 + biglen - 4, crc);
    PUT_32BIT(pkt->data + pktoffs, len);

    if (ssh->cipher)
	ssh->cipher->encrypt(ssh->v1_cipher_ctx,
			     pkt->data + pktoffs + 4, biglen);

    if (offset_p) *offset_p = pktoffs;
    return biglen + 4;		/* len(length+padding+type+data+CRC) */
}

static int s_write(Ssh ssh, void *data, int len)
{
    if (ssh->logctx)
	log_packet(ssh->logctx, PKT_OUTGOING, -1, NULL, data, len,
		   0, NULL, NULL);
    return sk_write(ssh->s, (char *)data, len);
}

static void s_wrpkt(Ssh ssh, struct Packet *pkt)
{
    int len, backlog, offset;
    len = s_wrpkt_prepare(ssh, pkt, &offset);
    backlog = s_write(ssh, pkt->data + offset, len);
    if (backlog > SSH_MAX_BACKLOG)
	ssh_throttle_all(ssh, 1, backlog);
    ssh_free_packet(pkt);
}

static void s_wrpkt_defer(Ssh ssh, struct Packet *pkt)
{
    int len, offset;
    len = s_wrpkt_prepare(ssh, pkt, &offset);
    if (ssh->deferred_len + len > ssh->deferred_size) {
	ssh->deferred_size = ssh->deferred_len + len + 128;
	ssh->deferred_send_data = sresize(ssh->deferred_send_data,
					  ssh->deferred_size,
					  unsigned char);
    }
    memcpy(ssh->deferred_send_data + ssh->deferred_len,
	   pkt->data + offset, len);
    ssh->deferred_len += len;
    ssh_free_packet(pkt);
}

/*
 * Construct a SSH-1 packet with the specified contents.
 * (This all-at-once interface used to be the only one, but now SSH-1
 * packets can also be constructed incrementally.)
 */
static struct Packet *construct_packet(Ssh ssh, int pkttype, va_list ap)
{
    int argtype;
    Bignum bn;
    struct Packet *pkt;

    pkt = ssh1_pkt_init(pkttype);

    while ((argtype = va_arg(ap, int)) != PKT_END) {
	unsigned char *argp, argchar;
	char *sargp;
	unsigned long argint;
	int arglen;
	switch (argtype) {
	  /* Actual fields in the packet */
	  case PKT_INT:
	    argint = va_arg(ap, int);
	    ssh_pkt_adduint32(pkt, argint);
	    break;
	  case PKT_CHAR:
	    argchar = (unsigned char) va_arg(ap, int);
	    ssh_pkt_addbyte(pkt, argchar);
	    break;
	  case PKT_DATA:
	    argp = va_arg(ap, unsigned char *);
	    arglen = va_arg(ap, int);
	    ssh_pkt_adddata(pkt, argp, arglen);
	    break;
	  case PKT_STR:
	    sargp = va_arg(ap, char *);
	    ssh_pkt_addstring(pkt, sargp);
	    break;
	  case PKT_BIGNUM:
	    bn = va_arg(ap, Bignum);
	    ssh1_pkt_addmp(pkt, bn);
	    break;
	  /* Tokens for modifications to packet logging */
	  case PKTT_PASSWORD:
	    dont_log_password(ssh, pkt, PKTLOG_BLANK);
	    break;
	  case PKTT_DATA:
	    dont_log_data(ssh, pkt, PKTLOG_OMIT);
	    break;
	  case PKTT_OTHER:
	    end_log_omission(ssh, pkt);
	    break;
	}
    }

    return pkt;
}

static void send_packet(Ssh ssh, int pkttype, ...)
{
    struct Packet *pkt;
    va_list ap;
    va_start(ap, pkttype);
    pkt = construct_packet(ssh, pkttype, ap);
    va_end(ap);
    s_wrpkt(ssh, pkt);
}

static void defer_packet(Ssh ssh, int pkttype, ...)
{
    struct Packet *pkt;
    va_list ap;
    va_start(ap, pkttype);
    pkt = construct_packet(ssh, pkttype, ap);
    va_end(ap);
    s_wrpkt_defer(ssh, pkt);
}

static int ssh_versioncmp(char *a, char *b)
{
    char *ae, *be;
    unsigned long av, bv;

    av = strtoul(a, &ae, 10);
    bv = strtoul(b, &be, 10);
    if (av != bv)
	return (av < bv ? -1 : +1);
    if (*ae == '.')
	ae++;
    if (*be == '.')
	be++;
    av = strtoul(ae, &ae, 10);
    bv = strtoul(be, &be, 10);
    if (av != bv)
	return (av < bv ? -1 : +1);
    return 0;
}

/*
 * Utility routines for putting an SSH-protocol `string' and
 * `uint32' into a hash state.
 */
static void hash_string(const struct ssh_hash *h, void *s, void *str, int len)
{
    unsigned char lenblk[4];
    PUT_32BIT(lenblk, len);
    h->bytes(s, lenblk, 4);
    h->bytes(s, str, len);
}

static void hash_uint32(const struct ssh_hash *h, void *s, unsigned i)
{
    unsigned char intblk[4];
    PUT_32BIT(intblk, i);
    h->bytes(s, intblk, 4);
}

/*
 * Packet construction functions. Mostly shared between SSH-1 and SSH-2.
 */
static void ssh_pkt_ensure(struct Packet *pkt, int length)
{
    if (pkt->maxlen < length) {
	unsigned char *body = pkt->body;
	int offset = body ? body - pkt->data : 0;
	pkt->maxlen = length + 256;
	pkt->data = sresize(pkt->data, pkt->maxlen + APIEXTRA, unsigned char);
	if (body) pkt->body = pkt->data + offset;
    }
}
static void ssh_pkt_adddata(struct Packet *pkt, void *data, int len)
{
    if (pkt->logmode != PKTLOG_EMIT) {
	pkt->nblanks++;
	pkt->blanks = sresize(pkt->blanks, pkt->nblanks, struct logblank_t);
	assert(pkt->body);
	pkt->blanks[pkt->nblanks-1].offset = pkt->length -
					     (pkt->body - pkt->data);
	pkt->blanks[pkt->nblanks-1].len = len;
	pkt->blanks[pkt->nblanks-1].type = pkt->logmode;
    }
    pkt->length += len;
    ssh_pkt_ensure(pkt, pkt->length);
    memcpy(pkt->data + pkt->length - len, data, len);
}
static void ssh_pkt_addbyte(struct Packet *pkt, unsigned char byte)
{
    ssh_pkt_adddata(pkt, &byte, 1);
}
static void ssh2_pkt_addbool(struct Packet *pkt, unsigned char value)
{
    ssh_pkt_adddata(pkt, &value, 1);
}
static void ssh_pkt_adduint32(struct Packet *pkt, unsigned long value)
{
    unsigned char x[4];
    PUT_32BIT(x, value);
    ssh_pkt_adddata(pkt, x, 4);
}
static void ssh_pkt_addstring_start(struct Packet *pkt)
{
    ssh_pkt_adduint32(pkt, 0);
    pkt->savedpos = pkt->length;
}
static void ssh_pkt_addstring_str(struct Packet *pkt, char *data)
{
    ssh_pkt_adddata(pkt, data, strlen(data));
    PUT_32BIT(pkt->data + pkt->savedpos - 4, pkt->length - pkt->savedpos);
}
static void ssh_pkt_addstring_data(struct Packet *pkt, char *data, int len)
{
    ssh_pkt_adddata(pkt, data, len);
    PUT_32BIT(pkt->data + pkt->savedpos - 4, pkt->length - pkt->savedpos);
}
static void ssh_pkt_addstring(struct Packet *pkt, char *data)
{
    ssh_pkt_addstring_start(pkt);
    ssh_pkt_addstring_str(pkt, data);
}
static void ssh1_pkt_addmp(struct Packet *pkt, Bignum b)
{
    int len = ssh1_bignum_length(b);
    unsigned char *data = snewn(len, unsigned char);
    (void) ssh1_write_bignum(data, b);
    ssh_pkt_adddata(pkt, data, len);
    sfree(data);
}
static unsigned char *ssh2_mpint_fmt(Bignum b, int *len)
{
    unsigned char *p;
    int i, n = (bignum_bitcount(b) + 7) / 8;
    p = snewn(n + 1, unsigned char);
    p[0] = 0;
    for (i = 1; i <= n; i++)
	p[i] = bignum_byte(b, n - i);
    i = 0;
    while (i <= n && p[i] == 0 && (p[i + 1] & 0x80) == 0)
	i++;
    memmove(p, p + i, n + 1 - i);
    *len = n + 1 - i;
    return p;
}
static void ssh2_pkt_addmp(struct Packet *pkt, Bignum b)
{
    unsigned char *p;
    int len;
    p = ssh2_mpint_fmt(b, &len);
    ssh_pkt_addstring_start(pkt);
    ssh_pkt_addstring_data(pkt, (char *)p, len);
    sfree(p);
}

static struct Packet *ssh1_pkt_init(int pkt_type)
{
    struct Packet *pkt = ssh_new_packet();
    pkt->length = 4 + 8;	    /* space for length + max padding */
    ssh_pkt_addbyte(pkt, pkt_type);
    pkt->body = pkt->data + pkt->length;
    return pkt;
}

/* For legacy code (SSH-1 and -2 packet construction used to be separate) */
#define ssh2_pkt_ensure(pkt, length) ssh_pkt_ensure(pkt, length)
#define ssh2_pkt_adddata(pkt, data, len) ssh_pkt_adddata(pkt, data, len)
#define ssh2_pkt_addbyte(pkt, byte) ssh_pkt_addbyte(pkt, byte)
#define ssh2_pkt_adduint32(pkt, value) ssh_pkt_adduint32(pkt, value)
#define ssh2_pkt_addstring_start(pkt) ssh_pkt_addstring_start(pkt)
#define ssh2_pkt_addstring_str(pkt, data) ssh_pkt_addstring_str(pkt, data)
#define ssh2_pkt_addstring_data(pkt, data, len) ssh_pkt_addstring_data(pkt, data, len)
#define ssh2_pkt_addstring(pkt, data) ssh_pkt_addstring(pkt, data)

static struct Packet *ssh2_pkt_init(int pkt_type)
{
    struct Packet *pkt = ssh_new_packet();
    pkt->length = 5; /* space for packet length + padding length */
    pkt->forcepad = 0;
    ssh_pkt_addbyte(pkt, (unsigned char) pkt_type);
    pkt->body = pkt->data + pkt->length; /* after packet type */
    return pkt;
}

/*
 * Construct an SSH-2 final-form packet: compress it, encrypt it,
 * put the MAC on it. Final packet, ready to be sent, is stored in
 * pkt->data. Total length is returned.
 */
static int ssh2_pkt_construct(Ssh ssh, struct Packet *pkt)
{
    int cipherblk, maclen, padding, i;

    if (ssh->logctx)
	log_packet(ssh->logctx, PKT_OUTGOING, pkt->data[5],
		   ssh2_pkt_type(ssh->pkt_kctx, ssh->pkt_actx, pkt->data[5]),
		   pkt->body, pkt->length - (pkt->body - pkt->data),
		   pkt->nblanks, pkt->blanks, &ssh->v2_outgoing_sequence);
    sfree(pkt->blanks); pkt->blanks = NULL;
    pkt->nblanks = 0;

    /*
     * Compress packet payload.
     */
    {
	unsigned char *newpayload;
	int newlen;
	if (ssh->cscomp &&
	    ssh->cscomp->compress(ssh->cs_comp_ctx, pkt->data + 5,
				  pkt->length - 5,
				  &newpayload, &newlen)) {
	    pkt->length = 5;
	    ssh2_pkt_adddata(pkt, newpayload, newlen);
	    sfree(newpayload);
	}
    }

    /*
     * Add padding. At least four bytes, and must also bring total
     * length (minus MAC) up to a multiple of the block size.
     * If pkt->forcepad is set, make sure the packet is at least that size
     * after padding.
     */
    cipherblk = ssh->cscipher ? ssh->cscipher->blksize : 8;  /* block size */
    cipherblk = cipherblk < 8 ? 8 : cipherblk;	/* or 8 if blksize < 8 */
    padding = 4;
    if (pkt->length + padding < pkt->forcepad)
	padding = pkt->forcepad - pkt->length;
    padding +=
	(cipherblk - (pkt->length + padding) % cipherblk) % cipherblk;
    assert(padding <= 255);
    maclen = ssh->csmac ? ssh->csmac->len : 0;
    ssh2_pkt_ensure(pkt, pkt->length + padding + maclen);
    pkt->data[4] = padding;
    for (i = 0; i < padding; i++)
	pkt->data[pkt->length + i] = random_byte();
    PUT_32BIT(pkt->data, pkt->length + padding - 4);
    if (ssh->csmac)
	ssh->csmac->generate(ssh->cs_mac_ctx, pkt->data,
			     pkt->length + padding,
			     ssh->v2_outgoing_sequence);
    ssh->v2_outgoing_sequence++;       /* whether or not we MACed */

    if (ssh->cscipher)
	ssh->cscipher->encrypt(ssh->cs_cipher_ctx,
			       pkt->data, pkt->length + padding);

    pkt->encrypted_len = pkt->length + padding;

    /* Ready-to-send packet starts at pkt->data. We return length. */
    return pkt->length + padding + maclen;
}

/*
 * Routines called from the main SSH code to send packets. There
 * are quite a few of these, because we have two separate
 * mechanisms for delaying the sending of packets:
 * 
 *  - In order to send an IGNORE message and a password message in
 *    a single fixed-length blob, we require the ability to
 *    concatenate the encrypted forms of those two packets _into_ a
 *    single blob and then pass it to our <network.h> transport
 *    layer in one go. Hence, there's a deferment mechanism which
 *    works after packet encryption.
 * 
 *  - In order to avoid sending any connection-layer messages
 *    during repeat key exchange, we have to queue up any such
 *    outgoing messages _before_ they are encrypted (and in
 *    particular before they're allocated sequence numbers), and
 *    then send them once we've finished.
 * 
 * I call these mechanisms `defer' and `queue' respectively, so as
 * to distinguish them reasonably easily.
 * 
 * The functions send_noqueue() and defer_noqueue() free the packet
 * structure they are passed. Every outgoing packet goes through
 * precisely one of these functions in its life; packets passed to
 * ssh2_pkt_send() or ssh2_pkt_defer() either go straight to one of
 * these or get queued, and then when the queue is later emptied
 * the packets are all passed to defer_noqueue().
 *
 * When using a CBC-mode cipher, it's necessary to ensure that an
 * attacker can't provide data to be encrypted using an IV that they
 * know.  We ensure this by prefixing each packet that might contain
 * user data with an SSH_MSG_IGNORE.  This is done using the deferral
 * mechanism, so in this case send_noqueue() ends up redirecting to
 * defer_noqueue().  If you don't like this inefficiency, don't use
 * CBC.
 */

static void ssh2_pkt_defer_noqueue(Ssh, struct Packet *, int);
static void ssh_pkt_defersend(Ssh);

/*
 * Send an SSH-2 packet immediately, without queuing or deferring.
 */
static void ssh2_pkt_send_noqueue(Ssh ssh, struct Packet *pkt)
{
    int len;
    int backlog;
    if (ssh->cscipher != NULL && (ssh->cscipher->flags & SSH_CIPHER_IS_CBC)) {
	/* We need to send two packets, so use the deferral mechanism. */
	ssh2_pkt_defer_noqueue(ssh, pkt, FALSE);
	ssh_pkt_defersend(ssh);
	return;
    }
    len = ssh2_pkt_construct(ssh, pkt);
    backlog = s_write(ssh, pkt->data, len);
    if (backlog > SSH_MAX_BACKLOG)
	ssh_throttle_all(ssh, 1, backlog);

    ssh->outgoing_data_size += pkt->encrypted_len;
    if (!ssh->kex_in_progress &&
	ssh->max_data_size != 0 &&
	ssh->outgoing_data_size > ssh->max_data_size)
	do_ssh2_transport(ssh, "too much data sent", -1, NULL);

    ssh_free_packet(pkt);
}

/*
 * Defer an SSH-2 packet.
 */
static void ssh2_pkt_defer_noqueue(Ssh ssh, struct Packet *pkt, int noignore)
{
    int len;
    if (ssh->cscipher != NULL && (ssh->cscipher->flags & SSH_CIPHER_IS_CBC) &&
	ssh->deferred_len == 0 && !noignore &&
	!(ssh->remote_bugs & BUG_CHOKES_ON_SSH2_IGNORE)) {
	/*
	 * Interpose an SSH_MSG_IGNORE to ensure that user data don't
	 * get encrypted with a known IV.
	 */
	struct Packet *ipkt = ssh2_pkt_init(SSH2_MSG_IGNORE);
	ssh2_pkt_addstring_start(ipkt);
	ssh2_pkt_defer_noqueue(ssh, ipkt, TRUE);
    }
    len = ssh2_pkt_construct(ssh, pkt);
    if (ssh->deferred_len + len > ssh->deferred_size) {
	ssh->deferred_size = ssh->deferred_len + len + 128;
	ssh->deferred_send_data = sresize(ssh->deferred_send_data,
					  ssh->deferred_size,
					  unsigned char);
    }
    memcpy(ssh->deferred_send_data + ssh->deferred_len, pkt->data, len);
    ssh->deferred_len += len;
    ssh->deferred_data_size += pkt->encrypted_len;
    ssh_free_packet(pkt);
}

/*
 * Queue an SSH-2 packet.
 */
static void ssh2_pkt_queue(Ssh ssh, struct Packet *pkt)
{
    assert(ssh->queueing);

    if (ssh->queuelen >= ssh->queuesize) {
	ssh->queuesize = ssh->queuelen + 32;
	ssh->queue = sresize(ssh->queue, ssh->queuesize, struct Packet *);
    }

    ssh->queue[ssh->queuelen++] = pkt;
}

/*
 * Either queue or send a packet, depending on whether queueing is
 * set.
 */
static void ssh2_pkt_send(Ssh ssh, struct Packet *pkt)
{
    if (ssh->queueing)
	ssh2_pkt_queue(ssh, pkt);
    else
	ssh2_pkt_send_noqueue(ssh, pkt);
}

/*
 * Either queue or defer a packet, depending on whether queueing is
 * set.
 */
static void ssh2_pkt_defer(Ssh ssh, struct Packet *pkt)
{
    if (ssh->queueing)
	ssh2_pkt_queue(ssh, pkt);
    else
	ssh2_pkt_defer_noqueue(ssh, pkt, FALSE);
}

/*
 * Send the whole deferred data block constructed by
 * ssh2_pkt_defer() or SSH-1's defer_packet().
 * 
 * The expected use of the defer mechanism is that you call
 * ssh2_pkt_defer() a few times, then call ssh_pkt_defersend(). If
 * not currently queueing, this simply sets up deferred_send_data
 * and then sends it. If we _are_ currently queueing, the calls to
 * ssh2_pkt_defer() put the deferred packets on to the queue
 * instead, and therefore ssh_pkt_defersend() has no deferred data
 * to send. Hence, there's no need to make it conditional on
 * ssh->queueing.
 */
static void ssh_pkt_defersend(Ssh ssh)
{
    int backlog;
    backlog = s_write(ssh, ssh->deferred_send_data, ssh->deferred_len);
    ssh->deferred_len = ssh->deferred_size = 0;
    sfree(ssh->deferred_send_data);
    ssh->deferred_send_data = NULL;
    if (backlog > SSH_MAX_BACKLOG)
	ssh_throttle_all(ssh, 1, backlog);

    ssh->outgoing_data_size += ssh->deferred_data_size;
    if (!ssh->kex_in_progress &&
	ssh->max_data_size != 0 &&
	ssh->outgoing_data_size > ssh->max_data_size)
	do_ssh2_transport(ssh, "too much data sent", -1, NULL);
    ssh->deferred_data_size = 0;
}

/*
 * Send a packet whose length needs to be disguised (typically
 * passwords or keyboard-interactive responses).
 */
static void ssh2_pkt_send_with_padding(Ssh ssh, struct Packet *pkt,
				       int padsize)
{
#if 0
    if (0) {
	/*
	 * The simplest way to do this is to adjust the
	 * variable-length padding field in the outgoing packet.
	 * 
	 * Currently compiled out, because some Cisco SSH servers
	 * don't like excessively padded packets (bah, why's it
	 * always Cisco?)
	 */
	pkt->forcepad = padsize;
	ssh2_pkt_send(ssh, pkt);
    } else
#endif
    {
	/*
	 * If we can't do that, however, an alternative approach is
	 * to use the pkt_defer mechanism to bundle the packet
	 * tightly together with an SSH_MSG_IGNORE such that their
	 * combined length is a constant. So first we construct the
	 * final form of this packet and defer its sending.
	 */
	ssh2_pkt_defer(ssh, pkt);

	/*
	 * Now construct an SSH_MSG_IGNORE which includes a string
	 * that's an exact multiple of the cipher block size. (If
	 * the cipher is NULL so that the block size is
	 * unavailable, we don't do this trick at all, because we
	 * gain nothing by it.)
	 */
	if (ssh->cscipher &&
	    !(ssh->remote_bugs & BUG_CHOKES_ON_SSH2_IGNORE)) {
	    int stringlen, i;

	    stringlen = (256 - ssh->deferred_len);
	    stringlen += ssh->cscipher->blksize - 1;
	    stringlen -= (stringlen % ssh->cscipher->blksize);
	    if (ssh->cscomp) {
		/*
		 * Temporarily disable actual compression, so we
		 * can guarantee to get this string exactly the
		 * length we want it. The compression-disabling
		 * routine should return an integer indicating how
		 * many bytes we should adjust our string length
		 * by.
		 */
		stringlen -=
		    ssh->cscomp->disable_compression(ssh->cs_comp_ctx);
	    }
	    pkt = ssh2_pkt_init(SSH2_MSG_IGNORE);
	    ssh2_pkt_addstring_start(pkt);
	    for (i = 0; i < stringlen; i++) {
		char c = (char) random_byte();
		ssh2_pkt_addstring_data(pkt, &c, 1);
	    }
	    ssh2_pkt_defer(ssh, pkt);
	}
	ssh_pkt_defersend(ssh);
    }
}

/*
 * Send all queued SSH-2 packets. We send them by means of
 * ssh2_pkt_defer_noqueue(), in case they included a pair of
 * packets that needed to be lumped together.
 */
static void ssh2_pkt_queuesend(Ssh ssh)
{
    int i;

    assert(!ssh->queueing);

    for (i = 0; i < ssh->queuelen; i++)
	ssh2_pkt_defer_noqueue(ssh, ssh->queue[i], FALSE);
    ssh->queuelen = 0;

    ssh_pkt_defersend(ssh);
}

#if 0
void bndebug(char *string, Bignum b)
{
    unsigned char *p;
    int i, len;
    p = ssh2_mpint_fmt(b, &len);
    debug(("%s", string));
    for (i = 0; i < len; i++)
	debug((" %02x", p[i]));
    debug(("\n"));
    sfree(p);
}
#endif

static void hash_mpint(const struct ssh_hash *h, void *s, Bignum b)
{
    unsigned char *p;
    int len;
    p = ssh2_mpint_fmt(b, &len);
    hash_string(h, s, p, len);
    sfree(p);
}

/*
 * Packet decode functions for both SSH-1 and SSH-2.
 */
static unsigned long ssh_pkt_getuint32(struct Packet *pkt)
{
    unsigned long value;
    if (pkt->length - pkt->savedpos < 4)
	return 0;		       /* arrgh, no way to decline (FIXME?) */
    value = GET_32BIT(pkt->body + pkt->savedpos);
    pkt->savedpos += 4;
    return value;
}
static int ssh2_pkt_getbool(struct Packet *pkt)
{
    unsigned long value;
    if (pkt->length - pkt->savedpos < 1)
	return 0;		       /* arrgh, no way to decline (FIXME?) */
    value = pkt->body[pkt->savedpos] != 0;
    pkt->savedpos++;
    return value;
}
static void ssh_pkt_getstring(struct Packet *pkt, char **p, int *length)
{
    int len;
    *p = NULL;
    *length = 0;
    if (pkt->length - pkt->savedpos < 4)
	return;
    len = GET_32BIT(pkt->body + pkt->savedpos);
    if (len < 0)
	return;
    *length = len;
    pkt->savedpos += 4;
    if (pkt->length - pkt->savedpos < *length)
	return;
    *p = (char *)(pkt->body + pkt->savedpos);
    pkt->savedpos += *length;
}
static void *ssh_pkt_getdata(struct Packet *pkt, int length)
{
    if (pkt->length - pkt->savedpos < length)
	return NULL;
    pkt->savedpos += length;
    return pkt->body + (pkt->savedpos - length);
}
static int ssh1_pkt_getrsakey(struct Packet *pkt, struct RSAKey *key,
			      unsigned char **keystr)
{
    int j;

    j = makekey(pkt->body + pkt->savedpos,
		pkt->length - pkt->savedpos,
		key, keystr, 0);

    if (j < 0)
	return FALSE;
    
    pkt->savedpos += j;
    assert(pkt->savedpos < pkt->length);

    return TRUE;
}
static Bignum ssh1_pkt_getmp(struct Packet *pkt)
{
    int j;
    Bignum b;

    j = ssh1_read_bignum(pkt->body + pkt->savedpos,
			 pkt->length - pkt->savedpos, &b);

    if (j < 0)
	return NULL;

    pkt->savedpos += j;
    return b;
}
static Bignum ssh2_pkt_getmp(struct Packet *pkt)
{
    char *p;
    int length;
    Bignum b;

    ssh_pkt_getstring(pkt, &p, &length);
    if (!p)
	return NULL;
    if (p[0] & 0x80)
	return NULL;
    b = bignum_from_bytes((unsigned char *)p, length);
    return b;
}

/*
 * Helper function to add an SSH-2 signature blob to a packet.
 * Expects to be shown the public key blob as well as the signature
 * blob. Normally works just like ssh2_pkt_addstring, but will
 * fiddle with the signature packet if necessary for
 * BUG_SSH2_RSA_PADDING.
 */
static void ssh2_add_sigblob(Ssh ssh, struct Packet *pkt,
			     void *pkblob_v, int pkblob_len,
			     void *sigblob_v, int sigblob_len)
{
    unsigned char *pkblob = (unsigned char *)pkblob_v;
    unsigned char *sigblob = (unsigned char *)sigblob_v;

    /* dmemdump(pkblob, pkblob_len); */
    /* dmemdump(sigblob, sigblob_len); */

    /*
     * See if this is in fact an ssh-rsa signature and a buggy
     * server; otherwise we can just do this the easy way.
     */
    if ((ssh->remote_bugs & BUG_SSH2_RSA_PADDING) &&
	(GET_32BIT(pkblob) == 7 && !memcmp(pkblob+4, "ssh-rsa", 7))) {
	int pos, len, siglen;

	/*
	 * Find the byte length of the modulus.
	 */

	pos = 4+7;		       /* skip over "ssh-rsa" */
	pos += 4 + GET_32BIT(pkblob+pos);   /* skip over exponent */
	len = GET_32BIT(pkblob+pos);   /* find length of modulus */
	pos += 4;		       /* find modulus itself */
	while (len > 0 && pkblob[pos] == 0)
	    len--, pos++;
	/* debug(("modulus length is %d\n", len)); */

	/*
	 * Now find the signature integer.
	 */
	pos = 4+7;		       /* skip over "ssh-rsa" */
	siglen = GET_32BIT(sigblob+pos);
	/* debug(("signature length is %d\n", siglen)); */

	if (len != siglen) {
	    unsigned char newlen[4];
	    ssh2_pkt_addstring_start(pkt);
	    ssh2_pkt_addstring_data(pkt, (char *)sigblob, pos);
	    /* dmemdump(sigblob, pos); */
	    pos += 4;		       /* point to start of actual sig */
	    PUT_32BIT(newlen, len);
	    ssh2_pkt_addstring_data(pkt, (char *)newlen, 4);
	    /* dmemdump(newlen, 4); */
	    newlen[0] = 0;
	    while (len-- > siglen) {
		ssh2_pkt_addstring_data(pkt, (char *)newlen, 1);
		/* dmemdump(newlen, 1); */
	    }
	    ssh2_pkt_addstring_data(pkt, (char *)(sigblob+pos), siglen);
	    /* dmemdump(sigblob+pos, siglen); */
	    return;
	}

	/* Otherwise fall through and do it the easy way. */
    }

    ssh2_pkt_addstring_start(pkt);
    ssh2_pkt_addstring_data(pkt, (char *)sigblob, sigblob_len);
}

/*
 * Examine the remote side's version string and compare it against
 * a list of known buggy implementations.
 */
static void ssh_detect_bugs(Ssh ssh, char *vstring)
{
    char *imp;			       /* pointer to implementation part */
    imp = vstring;
    imp += strcspn(imp, "-");
    if (*imp) imp++;
    imp += strcspn(imp, "-");
    if (*imp) imp++;

    ssh->remote_bugs = 0;

    /*
     * General notes on server version strings:
     *  - Not all servers reporting "Cisco-1.25" have all the bugs listed
     *    here -- in particular, we've heard of one that's perfectly happy
     *    with SSH1_MSG_IGNOREs -- but this string never seems to change,
     *    so we can't distinguish them.
     */
    if (ssh->cfg.sshbug_ignore1 == FORCE_ON ||
	(ssh->cfg.sshbug_ignore1 == AUTO &&
	 (!strcmp(imp, "1.2.18") || !strcmp(imp, "1.2.19") ||
	  !strcmp(imp, "1.2.20") || !strcmp(imp, "1.2.21") ||
	  !strcmp(imp, "1.2.22") || !strcmp(imp, "Cisco-1.25") ||
	  !strcmp(imp, "OSU_1.4alpha3") || !strcmp(imp, "OSU_1.5alpha4")))) {
	/*
	 * These versions don't support SSH1_MSG_IGNORE, so we have
	 * to use a different defence against password length
	 * sniffing.
	 */
	ssh->remote_bugs |= BUG_CHOKES_ON_SSH1_IGNORE;
	logevent("We believe remote version has SSH-1 ignore bug");
    }

    if (ssh->cfg.sshbug_plainpw1 == FORCE_ON ||
	(ssh->cfg.sshbug_plainpw1 == AUTO &&
	 (!strcmp(imp, "Cisco-1.25") || !strcmp(imp, "OSU_1.4alpha3")))) {
	/*
	 * These versions need a plain password sent; they can't
	 * handle having a null and a random length of data after
	 * the password.
	 */
	ssh->remote_bugs |= BUG_NEEDS_SSH1_PLAIN_PASSWORD;
	logevent("We believe remote version needs a plain SSH-1 password");
    }

    if (ssh->cfg.sshbug_rsa1 == FORCE_ON ||
	(ssh->cfg.sshbug_rsa1 == AUTO &&
	 (!strcmp(imp, "Cisco-1.25")))) {
	/*
	 * These versions apparently have no clue whatever about
	 * RSA authentication and will panic and die if they see
	 * an AUTH_RSA message.
	 */
	ssh->remote_bugs |= BUG_CHOKES_ON_RSA;
	logevent("We believe remote version can't handle SSH-1 RSA authentication");
    }

    if (ssh->cfg.sshbug_hmac2 == FORCE_ON ||
	(ssh->cfg.sshbug_hmac2 == AUTO &&
	 !wc_match("* VShell", imp) &&
	 (wc_match("2.1.0*", imp) || wc_match("2.0.*", imp) ||
	  wc_match("2.2.0*", imp) || wc_match("2.3.0*", imp) ||
	  wc_match("2.1 *", imp)))) {
	/*
	 * These versions have the HMAC bug.
	 */
	ssh->remote_bugs |= BUG_SSH2_HMAC;
	logevent("We believe remote version has SSH-2 HMAC bug");
    }

    if (ssh->cfg.sshbug_derivekey2 == FORCE_ON ||
	(ssh->cfg.sshbug_derivekey2 == AUTO &&
	 !wc_match("* VShell", imp) &&
	 (wc_match("2.0.0*", imp) || wc_match("2.0.10*", imp) ))) {
	/*
	 * These versions have the key-derivation bug (failing to
	 * include the literal shared secret in the hashes that
	 * generate the keys).
	 */
	ssh->remote_bugs |= BUG_SSH2_DERIVEKEY;
	logevent("We believe remote version has SSH-2 key-derivation bug");
    }

    if (ssh->cfg.sshbug_rsapad2 == FORCE_ON ||
	(ssh->cfg.sshbug_rsapad2 == AUTO &&
	 (wc_match("OpenSSH_2.[5-9]*", imp) ||
	  wc_match("OpenSSH_3.[0-2]*", imp)))) {
	/*
	 * These versions have the SSH-2 RSA padding bug.
	 */
	ssh->remote_bugs |= BUG_SSH2_RSA_PADDING;
	logevent("We believe remote version has SSH-2 RSA padding bug");
    }

    if (ssh->cfg.sshbug_pksessid2 == FORCE_ON ||
	(ssh->cfg.sshbug_pksessid2 == AUTO &&
	 wc_match("OpenSSH_2.[0-2]*", imp))) {
	/*
	 * These versions have the SSH-2 session-ID bug in
	 * public-key authentication.
	 */
	ssh->remote_bugs |= BUG_SSH2_PK_SESSIONID;
	logevent("We believe remote version has SSH-2 public-key-session-ID bug");
    }

    if (ssh->cfg.sshbug_rekey2 == FORCE_ON ||
	(ssh->cfg.sshbug_rekey2 == AUTO &&
	 (wc_match("DigiSSH_2.0", imp) ||
	  wc_match("OpenSSH_2.[0-4]*", imp) ||
	  wc_match("OpenSSH_2.5.[0-3]*", imp) ||
	  wc_match("Sun_SSH_1.0", imp) ||
	  wc_match("Sun_SSH_1.0.1", imp) ||
	  /* All versions <= 1.2.6 (they changed their format in 1.2.7) */
	  wc_match("WeOnlyDo-*", imp)))) {
	/*
	 * These versions have the SSH-2 rekey bug.
	 */
	ssh->remote_bugs |= BUG_SSH2_REKEY;
	logevent("We believe remote version has SSH-2 rekey bug");
    }

    if (ssh->cfg.sshbug_maxpkt2 == FORCE_ON ||
	(ssh->cfg.sshbug_maxpkt2 == AUTO &&
	 (wc_match("1.36_sshlib GlobalSCAPE", imp) ||
          wc_match("1.36 sshlib: GlobalScape", imp)))) {
	/*
	 * This version ignores our makpkt and needs to be throttled.
	 */
	ssh->remote_bugs |= BUG_SSH2_MAXPKT;
	logevent("We believe remote version ignores SSH-2 maximum packet size");
    }

    if (ssh->cfg.sshbug_ignore2 == FORCE_ON) {
	/*
	 * Servers that don't support SSH2_MSG_IGNORE. Currently,
	 * none detected automatically.
	 */
	ssh->remote_bugs |= BUG_CHOKES_ON_SSH2_IGNORE;
	logevent("We believe remote version has SSH-2 ignore bug");
    }
}

/*
 * The `software version' part of an SSH version string is required
 * to contain no spaces or minus signs.
 */
static void ssh_fix_verstring(char *str)
{
    /* Eat "SSH-<protoversion>-". */
    assert(*str == 'S'); str++;
    assert(*str == 'S'); str++;
    assert(*str == 'H'); str++;
    assert(*str == '-'); str++;
    while (*str && *str != '-') str++;
    assert(*str == '-'); str++;

    /* Convert minus signs and spaces in the remaining string into
     * underscores. */
    while (*str) {
        if (*str == '-' || *str == ' ')
            *str = '_';
        str++;
    }
}

/*
 * Send an appropriate SSH version string.
 */
static void ssh_send_verstring(Ssh ssh, char *svers)
{
    char *verstring;

    if (ssh->version == 2) {
	/*
	 * Construct a v2 version string.
	 */
	verstring = dupprintf("SSH-2.0-%s\015\012", sshver);
    } else {
	/*
	 * Construct a v1 version string.
	 */
	verstring = dupprintf("SSH-%s-%s\012",
			      (ssh_versioncmp(svers, "1.5") <= 0 ?
			       svers : "1.5"),
			      sshver);
    }

    ssh_fix_verstring(verstring);

    if (ssh->version == 2) {
	size_t len;
	/*
	 * Record our version string.
	 */
	len = strcspn(verstring, "\015\012");
	ssh->v_c = snewn(len + 1, char);
	memcpy(ssh->v_c, verstring, len);
	ssh->v_c[len] = 0;
    }

    logeventf(ssh, "We claim version: %.*s",
	      strcspn(verstring, "\015\012"), verstring);
    s_write(ssh, verstring, strlen(verstring));
    sfree(verstring);
}

static int do_ssh_init(Ssh ssh, unsigned char c)
{
    struct do_ssh_init_state {
	int vslen;
	char version[10];
	char *vstring;
	int vstrsize;
	int i;
	int proto1, proto2;
    };
    crState(do_ssh_init_state);

    crBegin(ssh->do_ssh_init_crstate);

    /* Search for a line beginning with the string "SSH-" in the input. */
    for (;;) {
	if (c != 'S') goto no;
	crReturn(1);
	if (c != 'S') goto no;
	crReturn(1);
	if (c != 'H') goto no;
	crReturn(1);
	if (c != '-') goto no;
	break;
      no:
	while (c != '\012')
	    crReturn(1);
	crReturn(1);
    }

    s->vstrsize = 16;
    s->vstring = snewn(s->vstrsize, char);
    strcpy(s->vstring, "SSH-");
    s->vslen = 4;
    s->i = 0;
    while (1) {
	crReturn(1);		       /* get another char */
	if (s->vslen >= s->vstrsize - 1) {
	    s->vstrsize += 16;
	    s->vstring = sresize(s->vstring, s->vstrsize, char);
	}
	s->vstring[s->vslen++] = c;
	if (s->i >= 0) {
	    if (c == '-') {
		s->version[s->i] = '\0';
		s->i = -1;
	    } else if (s->i < sizeof(s->version) - 1)
		s->version[s->i++] = c;
	} else if (c == '\012')
	    break;
    }

    ssh->agentfwd_enabled = FALSE;
    ssh->rdpkt2_state.incoming_sequence = 0;

    s->vstring[s->vslen] = 0;
    s->vstring[strcspn(s->vstring, "\015\012")] = '\0';/* remove EOL chars */
    logeventf(ssh, "Server version: %s", s->vstring);
    ssh_detect_bugs(ssh, s->vstring);

    /*
     * Decide which SSH protocol version to support.
     */

    /* Anything strictly below "2.0" means protocol 1 is supported. */
    s->proto1 = ssh_versioncmp(s->version, "2.0") < 0;
    /* Anything greater or equal to "1.99" means protocol 2 is supported. */
    s->proto2 = ssh_versioncmp(s->version, "1.99") >= 0;

    if (ssh->cfg.sshprot == 0 && !s->proto1) {
	bombout(("SSH protocol version 1 required by user but not provided by server"));
	crStop(0);
    }
    if (ssh->cfg.sshprot == 3 && !s->proto2) {
	bombout(("SSH protocol version 2 required by user but not provided by server"));
	crStop(0);
    }

    if (s->proto2 && (ssh->cfg.sshprot >= 2 || !s->proto1))
	ssh->version = 2;
    else
	ssh->version = 1;

    logeventf(ssh, "Using SSH protocol version %d", ssh->version);

    /* Send the version string, if we haven't already */
    if (ssh->cfg.sshprot != 3)
	ssh_send_verstring(ssh, s->version);

    if (ssh->version == 2) {
	size_t len;
	/*
	 * Record their version string.
	 */
	len = strcspn(s->vstring, "\015\012");
	ssh->v_s = snewn(len + 1, char);
	memcpy(ssh->v_s, s->vstring, len);
	ssh->v_s[len] = 0;
	    
	/*
	 * Initialise SSH-2 protocol.
	 */
	ssh->protocol = ssh2_protocol;
	ssh2_protocol_setup(ssh);
	ssh->s_rdpkt = ssh2_rdpkt;
    } else {
	/*
	 * Initialise SSH-1 protocol.
	 */
	ssh->protocol = ssh1_protocol;
	ssh1_protocol_setup(ssh);
	ssh->s_rdpkt = ssh1_rdpkt;
    }
    if (ssh->version == 2)
	do_ssh2_transport(ssh, NULL, -1, NULL);

    update_specials_menu(ssh->frontend);
    ssh->state = SSH_STATE_BEFORE_SIZE;
    ssh->pinger = pinger_new(&ssh->cfg, &ssh_backend, ssh);

    sfree(s->vstring);

    crFinish(0);
}

static void ssh_process_incoming_data(Ssh ssh,
				      unsigned char **data, int *datalen)
{
    struct Packet *pktin;

    pktin = ssh->s_rdpkt(ssh, data, datalen);
    if (pktin) {
	ssh->protocol(ssh, NULL, 0, pktin);
	ssh_free_packet(pktin);
    }
}

static void ssh_queue_incoming_data(Ssh ssh,
				    unsigned char **data, int *datalen)
{
    bufchain_add(&ssh->queued_incoming_data, *data, *datalen);
    *data += *datalen;
    *datalen = 0;
}

static void ssh_process_queued_incoming_data(Ssh ssh)
{
    void *vdata;
    unsigned char *data;
    int len, origlen;

    while (!ssh->frozen && bufchain_size(&ssh->queued_incoming_data)) {
	bufchain_prefix(&ssh->queued_incoming_data, &vdata, &len);
	data = vdata;
	origlen = len;

	while (!ssh->frozen && len > 0)
	    ssh_process_incoming_data(ssh, &data, &len);

	if (origlen > len)
	    bufchain_consume(&ssh->queued_incoming_data, origlen - len);
    }
}

static void ssh_set_frozen(Ssh ssh, int frozen)
{
    if (ssh->s)
	sk_set_frozen(ssh->s, frozen);
    ssh->frozen = frozen;
}

static void ssh_gotdata(Ssh ssh, unsigned char *data, int datalen)
{
    /* Log raw data, if we're in that mode. */
    if (ssh->logctx)
	log_packet(ssh->logctx, PKT_INCOMING, -1, NULL, data, datalen,
		   0, NULL, NULL);

    crBegin(ssh->ssh_gotdata_crstate);

    /*
     * To begin with, feed the characters one by one to the
     * protocol initialisation / selection function do_ssh_init().
     * When that returns 0, we're done with the initial greeting
     * exchange and can move on to packet discipline.
     */
    while (1) {
	int ret;		       /* need not be kept across crReturn */
	if (datalen == 0)
	    crReturnV;		       /* more data please */
	ret = do_ssh_init(ssh, *data);
	data++;
	datalen--;
	if (ret == 0)
	    break;
    }

    /*
     * We emerge from that loop when the initial negotiation is
     * over and we have selected an s_rdpkt function. Now pass
     * everything to s_rdpkt, and then pass the resulting packets
     * to the proper protocol handler.
     */

    while (1) {
	while (bufchain_size(&ssh->queued_incoming_data) > 0 || datalen > 0) {
	    if (ssh->frozen) {
		ssh_queue_incoming_data(ssh, &data, &datalen);
		/* This uses up all data and cannot cause anything interesting
		 * to happen; indeed, for anything to happen at all, we must
		 * return, so break out. */
		break;
	    } else if (bufchain_size(&ssh->queued_incoming_data) > 0) {
		/* This uses up some or all data, and may freeze the
		 * session. */
		ssh_process_queued_incoming_data(ssh);
	    } else {
		/* This uses up some or all data, and may freeze the
		 * session. */
		ssh_process_incoming_data(ssh, &data, &datalen);
	    }
	    /* FIXME this is probably EBW. */
	    if (ssh->state == SSH_STATE_CLOSED)
		return;
	}
	/* We're out of data. Go and get some more. */
	crReturnV;
    }
    crFinishV;
}

static int ssh_do_close(Ssh ssh, int notify_exit)
{
    int ret = 0;
    struct ssh_channel *c;

    ssh->state = SSH_STATE_CLOSED;
    expire_timer_context(ssh);
    if (ssh->s) {
        sk_close(ssh->s);
        ssh->s = NULL;
        if (notify_exit)
            notify_remote_exit(ssh->frontend);
        else
            ret = 1;
    }
    /*
     * Now we must shut down any port- and X-forwarded channels going
     * through this connection.
     */
    if (ssh->channels) {
	while (NULL != (c = index234(ssh->channels, 0))) {
	    switch (c->type) {
	      case CHAN_X11:
		x11_close(c->u.x11.s);
		break;
	      case CHAN_SOCKDATA:
	      case CHAN_SOCKDATA_DORMANT:
		pfd_close(c->u.pfd.s);
		break;
	    }
	    del234(ssh->channels, c); /* moving next one to index 0 */
	    if (ssh->version == 2)
		bufchain_clear(&c->v.v2.outbuffer);
	    sfree(c);
	}
    }
    /*
     * Go through port-forwardings, and close any associated
     * listening sockets.
     */
    if (ssh->portfwds) {
	struct ssh_portfwd *pf;
	while (NULL != (pf = index234(ssh->portfwds, 0))) {
	    /* Dispose of any listening socket. */
	    if (pf->local)
		pfd_terminate(pf->local);
	    del234(ssh->portfwds, pf); /* moving next one to index 0 */
	    free_portfwd(pf);
	}
	freetree234(ssh->portfwds);
	ssh->portfwds = NULL;
    }

    return ret;
}

static void ssh_log(Plug plug, int type, SockAddr addr, int port,
		    const char *error_msg, int error_code)
{
    Ssh ssh = (Ssh) plug;
    char addrbuf[256], *msg;

    sk_getaddr(addr, addrbuf, lenof(addrbuf));

    if (type == 0)
	msg = dupprintf("Connecting to %s port %d", addrbuf, port);
    else
	msg = dupprintf("Failed to connect to %s: %s", addrbuf, error_msg);

    logevent(msg);
    sfree(msg);
}

static int ssh_closing(Plug plug, const char *error_msg, int error_code,
		       int calling_back)
{
    Ssh ssh = (Ssh) plug;
    int need_notify = ssh_do_close(ssh, FALSE);

    if (!error_msg) {
	if (!ssh->close_expected)
	    error_msg = "Server unexpectedly closed network connection";
	else
	    error_msg = "Server closed network connection";
    }

    if (ssh->close_expected && ssh->clean_exit && ssh->exitcode < 0)
	ssh->exitcode = 0;

    if (need_notify)
        notify_remote_exit(ssh->frontend);

    if (error_msg)
	logevent(error_msg);
    if (!ssh->close_expected || !ssh->clean_exit)
	connection_fatal(ssh->frontend, "%s", error_msg);
    return 0;
}

static int ssh_receive(Plug plug, int urgent, char *data, int len)
{
    Ssh ssh = (Ssh) plug;
    ssh_gotdata(ssh, (unsigned char *)data, len);
    if (ssh->state == SSH_STATE_CLOSED) {
	ssh_do_close(ssh, TRUE);
	return 0;
    }
    return 1;
}

static void ssh_sent(Plug plug, int bufsize)
{
    Ssh ssh = (Ssh) plug;
    /*
     * If the send backlog on the SSH socket itself clears, we
     * should unthrottle the whole world if it was throttled.
     */
    if (bufsize < SSH_MAX_BACKLOG)
	ssh_throttle_all(ssh, 0, bufsize);
}

/*
 * Connect to specified host and port.
 * Returns an error message, or NULL on success.
 * Also places the canonical host name into `realhost'. It must be
 * freed by the caller.
 */
static const char *connect_to_host(Ssh ssh, char *host, int port,
				   char **realhost, int nodelay, int keepalive)
{
    static const struct plug_function_table fn_table = {
	ssh_log,
	ssh_closing,
	ssh_receive,
	ssh_sent,
	NULL
    };

    SockAddr addr;
    const char *err;

    if (*ssh->cfg.loghost) {
	char *colon;

	ssh->savedhost = dupstr(ssh->cfg.loghost);
	ssh->savedport = 22;	       /* default ssh port */

	/*
	 * A colon suffix on savedhost also lets us affect
	 * savedport.
	 * 
	 * (FIXME: do something about IPv6 address literals here.)
	 */
	colon = strrchr(ssh->savedhost, ':');
	if (colon) {
	    *colon++ = '\0';
	    if (*colon)
		ssh->savedport = atoi(colon);
	}
    } else {
	ssh->savedhost = dupstr(host);
	if (port < 0)
	    port = 22;		       /* default ssh port */
	ssh->savedport = port;
    }

    /*
     * Try to find host.
     */
    logeventf(ssh, "Looking up host \"%s\"%s", host,
	      (ssh->cfg.addressfamily == ADDRTYPE_IPV4 ? " (IPv4)" :
	       (ssh->cfg.addressfamily == ADDRTYPE_IPV6 ? " (IPv6)" : "")));
    addr = name_lookup(host, port, realhost, &ssh->cfg,
		       ssh->cfg.addressfamily);
    if ((err = sk_addr_error(addr)) != NULL) {
	sk_addr_free(addr);
	return err;
    }
    ssh->fullhostname = dupstr(*realhost);   /* save in case of GSSAPI */

    /*
     * Open socket.
     */
    ssh->fn = &fn_table;
    ssh->s = new_connection(addr, *realhost, port,
			    0, 1, nodelay, keepalive, (Plug) ssh, &ssh->cfg);
    if ((err = sk_socket_error(ssh->s)) != NULL) {
	ssh->s = NULL;
	notify_remote_exit(ssh->frontend);
	return err;
    }

    /*
     * If the SSH version number's fixed, set it now, and if it's SSH-2,
     * send the version string too.
     */
    if (ssh->cfg.sshprot == 0)
	ssh->version = 1;
    if (ssh->cfg.sshprot == 3) {
	ssh->version = 2;
	ssh_send_verstring(ssh, NULL);
    }

    /*
     * loghost, if configured, overrides realhost.
     */
    if (*ssh->cfg.loghost) {
	sfree(*realhost);
	*realhost = dupstr(ssh->cfg.loghost);
    }

    return NULL;
}

/*
 * Throttle or unthrottle the SSH connection.
 */
static void ssh_throttle_conn(Ssh ssh, int adjust)
{
    int old_count = ssh->conn_throttle_count;
    ssh->conn_throttle_count += adjust;
    assert(ssh->conn_throttle_count >= 0);
    if (ssh->conn_throttle_count && !old_count) {
	ssh_set_frozen(ssh, 1);
    } else if (!ssh->conn_throttle_count && old_count) {
	ssh_set_frozen(ssh, 0);
    }
}

/*
 * Throttle or unthrottle _all_ local data streams (for when sends
 * on the SSH connection itself back up).
 */
static void ssh_throttle_all(Ssh ssh, int enable, int bufsize)
{
    int i;
    struct ssh_channel *c;

    if (enable == ssh->throttled_all)
	return;
    ssh->throttled_all = enable;
    ssh->overall_bufsize = bufsize;
    if (!ssh->channels)
	return;
    for (i = 0; NULL != (c = index234(ssh->channels, i)); i++) {
	switch (c->type) {
	  case CHAN_MAINSESSION:
	    /*
	     * This is treated separately, outside the switch.
	     */
	    break;
	  case CHAN_X11:
	    x11_override_throttle(c->u.x11.s, enable);
	    break;
	  case CHAN_AGENT:
	    /* Agent channels require no buffer management. */
	    break;
	  case CHAN_SOCKDATA:
	    pfd_override_throttle(c->u.pfd.s, enable);
	    break;
	}
    }
}

static void ssh_agent_callback(void *sshv, void *reply, int replylen)
{
    Ssh ssh = (Ssh) sshv;

    ssh->agent_response = reply;
    ssh->agent_response_len = replylen;

    if (ssh->version == 1)
	do_ssh1_login(ssh, NULL, -1, NULL);
    else
	do_ssh2_authconn(ssh, NULL, -1, NULL);
}

static void ssh_dialog_callback(void *sshv, int ret)
{
    Ssh ssh = (Ssh) sshv;

    ssh->user_response = ret;

    if (ssh->version == 1)
	do_ssh1_login(ssh, NULL, -1, NULL);
    else
	do_ssh2_transport(ssh, NULL, -1, NULL);

    /*
     * This may have unfrozen the SSH connection, so do a
     * queued-data run.
     */
    ssh_process_queued_incoming_data(ssh);
}

static void ssh_agentf_callback(void *cv, void *reply, int replylen)
{
    struct ssh_channel *c = (struct ssh_channel *)cv;
    Ssh ssh = c->ssh;
    void *sentreply = reply;

    if (!sentreply) {
	/* Fake SSH_AGENT_FAILURE. */
	sentreply = "\0\0\0\1\5";
	replylen = 5;
    }
    if (ssh->version == 2) {
	ssh2_add_channel_data(c, sentreply, replylen);
	ssh2_try_send(c);
    } else {
	send_packet(ssh, SSH1_MSG_CHANNEL_DATA,
		    PKT_INT, c->remoteid,
		    PKT_INT, replylen,
		    PKTT_DATA,
		    PKT_DATA, sentreply, replylen,
		    PKTT_OTHER,
		    PKT_END);
    }
    if (reply)
	sfree(reply);
}

/*
 * Client-initiated disconnection. Send a DISCONNECT if `wire_reason'
 * non-NULL, otherwise just close the connection. `client_reason' == NULL
 * => log `wire_reason'.
 */
static void ssh_disconnect(Ssh ssh, char *client_reason, char *wire_reason,
			   int code, int clean_exit)
{
    char *error;
    if (!client_reason)
	client_reason = wire_reason;
    if (client_reason)
	error = dupprintf("Disconnected: %s", client_reason);
    else
	error = dupstr("Disconnected");
    if (wire_reason) {
	if (ssh->version == 1) {
	    send_packet(ssh, SSH1_MSG_DISCONNECT, PKT_STR, wire_reason,
			PKT_END);
	} else if (ssh->version == 2) {
	    struct Packet *pktout = ssh2_pkt_init(SSH2_MSG_DISCONNECT);
	    ssh2_pkt_adduint32(pktout, code);
	    ssh2_pkt_addstring(pktout, wire_reason);
	    ssh2_pkt_addstring(pktout, "en");	/* language tag */
	    ssh2_pkt_send_noqueue(ssh, pktout);
	}
    }
    ssh->close_expected = TRUE;
    ssh->clean_exit = clean_exit;
    ssh_closing((Plug)ssh, error, 0, 0);
    sfree(error);
}

/*
 * Handle the key exchange and user authentication phases.
 */
static int do_ssh1_login(Ssh ssh, unsigned char *in, int inlen,
			 struct Packet *pktin)
{
    int i, j, ret;
    unsigned char cookie[8], *ptr;
    struct RSAKey servkey, hostkey;
    struct MD5Context md5c;
    struct do_ssh1_login_state {
	int len;
	unsigned char *rsabuf, *keystr1, *keystr2;
	unsigned long supported_ciphers_mask, supported_auths_mask;
	int tried_publickey, tried_agent;
	int tis_auth_refused, ccard_auth_refused;
	unsigned char session_id[16];
	int cipher_type;
	char username[100];
	void *publickey_blob;
	int publickey_bloblen;
	char *publickey_comment;
	int publickey_encrypted;
	prompts_t *cur_prompt;
	char c;
	int pwpkt_type;
	unsigned char request[5], *response, *p;
	int responselen;
	int keyi, nkeys;
	int authed;
	struct RSAKey key;
	Bignum challenge;
	char *commentp;
	int commentlen;
        int dlgret;
    };
    crState(do_ssh1_login_state);

    crBegin(ssh->do_ssh1_login_crstate);

    if (!pktin)
	crWaitUntil(pktin);

    if (pktin->type != SSH1_SMSG_PUBLIC_KEY) {
	bombout(("Public key packet not received"));
	crStop(0);
    }

    logevent("Received public keys");

    ptr = ssh_pkt_getdata(pktin, 8);
    if (!ptr) {
	bombout(("SSH-1 public key packet stopped before random cookie"));
	crStop(0);
    }
    memcpy(cookie, ptr, 8);

    if (!ssh1_pkt_getrsakey(pktin, &servkey, &s->keystr1) ||
	!ssh1_pkt_getrsakey(pktin, &hostkey, &s->keystr2)) {	
	bombout(("Failed to read SSH-1 public keys from public key packet"));
	crStop(0);
    }

    /*
     * Log the host key fingerprint.
     */
    {
	char logmsg[80];
	logevent("Host key fingerprint is:");
	strcpy(logmsg, "      ");
	hostkey.comment = NULL;
	rsa_fingerprint(logmsg + strlen(logmsg),
			sizeof(logmsg) - strlen(logmsg), &hostkey);
	logevent(logmsg);
    }

    ssh->v1_remote_protoflags = ssh_pkt_getuint32(pktin);
    s->supported_ciphers_mask = ssh_pkt_getuint32(pktin);
    s->supported_auths_mask = ssh_pkt_getuint32(pktin);
    if ((ssh->remote_bugs & BUG_CHOKES_ON_RSA))
	s->supported_auths_mask &= ~(1 << SSH1_AUTH_RSA);

    ssh->v1_local_protoflags =
	ssh->v1_remote_protoflags & SSH1_PROTOFLAGS_SUPPORTED;
    ssh->v1_local_protoflags |= SSH1_PROTOFLAG_SCREEN_NUMBER;

    MD5Init(&md5c);
    MD5Update(&md5c, s->keystr2, hostkey.bytes);
    MD5Update(&md5c, s->keystr1, servkey.bytes);
    MD5Update(&md5c, cookie, 8);
    MD5Final(s->session_id, &md5c);

    for (i = 0; i < 32; i++)
	ssh->session_key[i] = random_byte();

    /*
     * Verify that the `bits' and `bytes' parameters match.
     */
    if (hostkey.bits > hostkey.bytes * 8 ||
	servkey.bits > servkey.bytes * 8) {
	bombout(("SSH-1 public keys were badly formatted"));
	crStop(0);
    }

    s->len = (hostkey.bytes > servkey.bytes ? hostkey.bytes : servkey.bytes);

    s->rsabuf = snewn(s->len, unsigned char);

    /*
     * Verify the host key.
     */
    {
	/*
	 * First format the key into a string.
	 */
	int len = rsastr_len(&hostkey);
	char fingerprint[100];
	char *keystr = snewn(len, char);
	rsastr_fmt(keystr, &hostkey);
	rsa_fingerprint(fingerprint, sizeof(fingerprint), &hostkey);

        ssh_set_frozen(ssh, 1);
	s->dlgret = verify_ssh_host_key(ssh->frontend,
                                        ssh->savedhost, ssh->savedport,
                                        "rsa", keystr, fingerprint,
                                        ssh_dialog_callback, ssh);
	sfree(keystr);
        if (s->dlgret < 0) {
            do {
                crReturn(0);
                if (pktin) {
                    bombout(("Unexpected data from server while waiting"
                             " for user host key response"));
                    crStop(0);
                }
            } while (pktin || inlen > 0);
            s->dlgret = ssh->user_response;
        }
        ssh_set_frozen(ssh, 0);

        if (s->dlgret == 0) {
	    ssh_disconnect(ssh, "User aborted at host key verification",
			   NULL, 0, TRUE);
	    crStop(0);
        }
    }

    for (i = 0; i < 32; i++) {
	s->rsabuf[i] = ssh->session_key[i];
	if (i < 16)
	    s->rsabuf[i] ^= s->session_id[i];
    }

    if (hostkey.bytes > servkey.bytes) {
	ret = rsaencrypt(s->rsabuf, 32, &servkey);
	if (ret)
	    ret = rsaencrypt(s->rsabuf, servkey.bytes, &hostkey);
    } else {
	ret = rsaencrypt(s->rsabuf, 32, &hostkey);
	if (ret)
	    ret = rsaencrypt(s->rsabuf, hostkey.bytes, &servkey);
    }
    if (!ret) {
	bombout(("SSH-1 public key encryptions failed due to bad formatting"));
	crStop(0);	
    }

    logevent("Encrypted session key");

    {
	int cipher_chosen = 0, warn = 0;
	char *cipher_string = NULL;
	int i;
	for (i = 0; !cipher_chosen && i < CIPHER_MAX; i++) {
	    int next_cipher = ssh->cfg.ssh_cipherlist[i];
	    if (next_cipher == CIPHER_WARN) {
		/* If/when we choose a cipher, warn about it */
		warn = 1;
	    } else if (next_cipher == CIPHER_AES) {
		/* XXX Probably don't need to mention this. */
		logevent("AES not supported in SSH-1, skipping");
	    } else {
		switch (next_cipher) {
		  case CIPHER_3DES:     s->cipher_type = SSH_CIPHER_3DES;
					cipher_string = "3DES"; break;
		  case CIPHER_BLOWFISH: s->cipher_type = SSH_CIPHER_BLOWFISH;
					cipher_string = "Blowfish"; break;
		  case CIPHER_DES:	s->cipher_type = SSH_CIPHER_DES;
					cipher_string = "single-DES"; break;
		}
		if (s->supported_ciphers_mask & (1 << s->cipher_type))
		    cipher_chosen = 1;
	    }
	}
	if (!cipher_chosen) {
	    if ((s->supported_ciphers_mask & (1 << SSH_CIPHER_3DES)) == 0)
		bombout(("Server violates SSH-1 protocol by not "
			 "supporting 3DES encryption"));
	    else
		/* shouldn't happen */
		bombout(("No supported ciphers found"));
	    crStop(0);
	}

	/* Warn about chosen cipher if necessary. */
	if (warn) {
            ssh_set_frozen(ssh, 1);
	    s->dlgret = askalg(ssh->frontend, "cipher", cipher_string,
			       ssh_dialog_callback, ssh);
	    if (s->dlgret < 0) {
		do {
		    crReturn(0);
		    if (pktin) {
			bombout(("Unexpected data from server while waiting"
				 " for user response"));
			crStop(0);
		    }
		} while (pktin || inlen > 0);
		s->dlgret = ssh->user_response;
	    }
            ssh_set_frozen(ssh, 0);
	    if (s->dlgret == 0) {
		ssh_disconnect(ssh, "User aborted at cipher warning", NULL,
			       0, TRUE);
		crStop(0);
	    }
        }
    }

    switch (s->cipher_type) {
      case SSH_CIPHER_3DES:
	logevent("Using 3DES encryption");
	break;
      case SSH_CIPHER_DES:
	logevent("Using single-DES encryption");
	break;
      case SSH_CIPHER_BLOWFISH:
	logevent("Using Blowfish encryption");
	break;
    }

    send_packet(ssh, SSH1_CMSG_SESSION_KEY,
		PKT_CHAR, s->cipher_type,
		PKT_DATA, cookie, 8,
		PKT_CHAR, (s->len * 8) >> 8, PKT_CHAR, (s->len * 8) & 0xFF,
		PKT_DATA, s->rsabuf, s->len,
		PKT_INT, ssh->v1_local_protoflags, PKT_END);

    logevent("Trying to enable encryption...");

    sfree(s->rsabuf);

    ssh->cipher = (s->cipher_type == SSH_CIPHER_BLOWFISH ? &ssh_blowfish_ssh1 :
		   s->cipher_type == SSH_CIPHER_DES ? &ssh_des :
		   &ssh_3des);
    ssh->v1_cipher_ctx = ssh->cipher->make_context();
    ssh->cipher->sesskey(ssh->v1_cipher_ctx, ssh->session_key);
    logeventf(ssh, "Initialised %s encryption", ssh->cipher->text_name);

    ssh->crcda_ctx = crcda_make_context();
    logevent("Installing CRC compensation attack detector");

    if (servkey.modulus) {
	sfree(servkey.modulus);
	servkey.modulus = NULL;
    }
    if (servkey.exponent) {
	sfree(servkey.exponent);
	servkey.exponent = NULL;
    }
    if (hostkey.modulus) {
	sfree(hostkey.modulus);
	hostkey.modulus = NULL;
    }
    if (hostkey.exponent) {
	sfree(hostkey.exponent);
	hostkey.exponent = NULL;
    }
    crWaitUntil(pktin);

    if (pktin->type != SSH1_SMSG_SUCCESS) {
	bombout(("Encryption not successfully enabled"));
	crStop(0);
    }

    logevent("Successfully started encryption");

    fflush(stdout); /* FIXME eh? */
    {
	if (!get_remote_username(&ssh->cfg, s->username,
				 sizeof(s->username))) {
	    int ret; /* need not be kept over crReturn */
	    s->cur_prompt = new_prompts(ssh->frontend);
	    s->cur_prompt->to_server = TRUE;
	    s->cur_prompt->name = dupstr("SSH login name");
	    add_prompt(s->cur_prompt, dupstr("login as: "), TRUE,
		       lenof(s->username)); 
	    ret = get_userpass_input(s->cur_prompt, NULL, 0);
	    while (ret < 0) {
		ssh->send_ok = 1;
		crWaitUntil(!pktin);
		ret = get_userpass_input(s->cur_prompt, in, inlen);
		ssh->send_ok = 0;
	    }
	    if (!ret) {
		/*
		 * Failed to get a username. Terminate.
		 */
		free_prompts(s->cur_prompt);
		ssh_disconnect(ssh, "No username provided", NULL, 0, TRUE);
		crStop(0);
	    }
	    memcpy(s->username, s->cur_prompt->prompts[0]->result,
		   lenof(s->username));
	    free_prompts(s->cur_prompt);
	}

	send_packet(ssh, SSH1_CMSG_USER, PKT_STR, s->username, PKT_END);
	{
	    char *userlog = dupprintf("Sent username \"%s\"", s->username);
	    logevent(userlog);
	    if (flags & FLAG_INTERACTIVE &&
		(!((flags & FLAG_STDERR) && (flags & FLAG_VERBOSE)))) {
		c_write_str(ssh, userlog);
		c_write_str(ssh, "\r\n");
	    }
	    sfree(userlog);
	}
    }

    crWaitUntil(pktin);

    if ((s->supported_auths_mask & (1 << SSH1_AUTH_RSA)) == 0) {
	/* We must not attempt PK auth. Pretend we've already tried it. */
	s->tried_publickey = s->tried_agent = 1;
    } else {
	s->tried_publickey = s->tried_agent = 0;
    }
    s->tis_auth_refused = s->ccard_auth_refused = 0;
    /*
     * Load the public half of any configured keyfile for later use.
     */
    if (!filename_is_null(ssh->cfg.keyfile)) {
	int keytype;
	logeventf(ssh, "Reading private key file \"%.150s\"",
		  filename_to_str(&ssh->cfg.keyfile));
	keytype = key_type(&ssh->cfg.keyfile);
	if (keytype == SSH_KEYTYPE_SSH1) {
	    const char *error;
	    if (rsakey_pubblob(&ssh->cfg.keyfile,
			       &s->publickey_blob, &s->publickey_bloblen,
			       &s->publickey_comment, &error)) {
		s->publickey_encrypted = rsakey_encrypted(&ssh->cfg.keyfile,
							  NULL);
	    } else {
		char *msgbuf;
		logeventf(ssh, "Unable to load private key (%s)", error);
		msgbuf = dupprintf("Unable to load private key file "
				   "\"%.150s\" (%s)\r\n",
				   filename_to_str(&ssh->cfg.keyfile),
				   error);
		c_write_str(ssh, msgbuf);
		sfree(msgbuf);
		s->publickey_blob = NULL;
	    }
	} else {
	    char *msgbuf;
	    logeventf(ssh, "Unable to use this key file (%s)",
		      key_type_to_str(keytype));
	    msgbuf = dupprintf("Unable to use key file \"%.150s\""
			       " (%s)\r\n",
			       filename_to_str(&ssh->cfg.keyfile),
			       key_type_to_str(keytype));
	    c_write_str(ssh, msgbuf);
	    sfree(msgbuf);
	    s->publickey_blob = NULL;
	}
    } else
	s->publickey_blob = NULL;

    while (pktin->type == SSH1_SMSG_FAILURE) {
	s->pwpkt_type = SSH1_CMSG_AUTH_PASSWORD;

	if (ssh->cfg.tryagent && agent_exists() && !s->tried_agent) {
	    /*
	     * Attempt RSA authentication using Pageant.
	     */
	    void *r;

	    s->authed = FALSE;
	    s->tried_agent = 1;
	    logevent("Pageant is running. Requesting keys.");

	    /* Request the keys held by the agent. */
	    PUT_32BIT(s->request, 1);
	    s->request[4] = SSH1_AGENTC_REQUEST_RSA_IDENTITIES;
	    if (!agent_query(s->request, 5, &r, &s->responselen,
			     ssh_agent_callback, ssh)) {
		do {
		    crReturn(0);
		    if (pktin) {
			bombout(("Unexpected data from server while waiting"
				 " for agent response"));
			crStop(0);
		    }
		} while (pktin || inlen > 0);
		r = ssh->agent_response;
		s->responselen = ssh->agent_response_len;
	    }
	    s->response = (unsigned char *) r;
	    if (s->response && s->responselen >= 5 &&
		s->response[4] == SSH1_AGENT_RSA_IDENTITIES_ANSWER) {
		s->p = s->response + 5;
		s->nkeys = GET_32BIT(s->p);
		s->p += 4;
		logeventf(ssh, "Pageant has %d SSH-1 keys", s->nkeys);
		for (s->keyi = 0; s->keyi < s->nkeys; s->keyi++) {
		    unsigned char *pkblob = s->p;
		    s->p += 4;
		    {
			int n, ok = FALSE;
			do {	       /* do while (0) to make breaking easy */
			    n = ssh1_read_bignum
				(s->p, s->responselen-(s->p-s->response),
				 &s->key.exponent);
			    if (n < 0)
				break;
			    s->p += n;
			    n = ssh1_read_bignum
				(s->p, s->responselen-(s->p-s->response),
				 &s->key.modulus);
			    if (n < 0)
			    break;
			    s->p += n;
			    if (s->responselen - (s->p-s->response) < 4)
				break;
			    s->commentlen = GET_32BIT(s->p);
			    s->p += 4;
			    if (s->responselen - (s->p-s->response) <
				s->commentlen)
				break;
			    s->commentp = (char *)s->p;
			    s->p += s->commentlen;
			    ok = TRUE;
			} while (0);
			if (!ok) {
			    logevent("Pageant key list packet was truncated");
			    break;
			}
		    }
		    if (s->publickey_blob) {
			if (!memcmp(pkblob, s->publickey_blob,
				    s->publickey_bloblen)) {
			    logeventf(ssh, "Pageant key #%d matches "
				      "configured key file", s->keyi);
			    s->tried_publickey = 1;
			} else
			    /* Skip non-configured key */
			    continue;
		    }
		    logeventf(ssh, "Trying Pageant key #%d", s->keyi);
		    send_packet(ssh, SSH1_CMSG_AUTH_RSA,
				PKT_BIGNUM, s->key.modulus, PKT_END);
		    crWaitUntil(pktin);
		    if (pktin->type != SSH1_SMSG_AUTH_RSA_CHALLENGE) {
			logevent("Key refused");
			continue;
		    }
		    logevent("Received RSA challenge");
		    if ((s->challenge = ssh1_pkt_getmp(pktin)) == NULL) {
			bombout(("Server's RSA challenge was badly formatted"));
			crStop(0);
		    }

		    {
			char *agentreq, *q, *ret;
			void *vret;
			int len, retlen;
			len = 1 + 4;   /* message type, bit count */
			len += ssh1_bignum_length(s->key.exponent);
			len += ssh1_bignum_length(s->key.modulus);
			len += ssh1_bignum_length(s->challenge);
			len += 16;     /* session id */
			len += 4;      /* response format */
			agentreq = snewn(4 + len, char);
			PUT_32BIT(agentreq, len);
			q = agentreq + 4;
			*q++ = SSH1_AGENTC_RSA_CHALLENGE;
			PUT_32BIT(q, bignum_bitcount(s->key.modulus));
			q += 4;
			q += ssh1_write_bignum(q, s->key.exponent);
			q += ssh1_write_bignum(q, s->key.modulus);
			q += ssh1_write_bignum(q, s->challenge);
			memcpy(q, s->session_id, 16);
			q += 16;
			PUT_32BIT(q, 1);	/* response format */
			if (!agent_query(agentreq, len + 4, &vret, &retlen,
					 ssh_agent_callback, ssh)) {
			    sfree(agentreq);
			    do {
				crReturn(0);
				if (pktin) {
				    bombout(("Unexpected data from server"
					     " while waiting for agent"
					     " response"));
				    crStop(0);
				}
			    } while (pktin || inlen > 0);
			    vret = ssh->agent_response;
			    retlen = ssh->agent_response_len;
			} else
			    sfree(agentreq);
			ret = vret;
			if (ret) {
			    if (ret[4] == SSH1_AGENT_RSA_RESPONSE) {
				logevent("Sending Pageant's response");
				send_packet(ssh, SSH1_CMSG_AUTH_RSA_RESPONSE,
					    PKT_DATA, ret + 5, 16,
					    PKT_END);
				sfree(ret);
				crWaitUntil(pktin);
				if (pktin->type == SSH1_SMSG_SUCCESS) {
				    logevent
					("Pageant's response accepted");
				    if (flags & FLAG_VERBOSE) {
					c_write_str(ssh, "Authenticated using"
						    " RSA key \"");
					c_write(ssh, s->commentp,
						s->commentlen);
					c_write_str(ssh, "\" from agent\r\n");
				    }
				    s->authed = TRUE;
				} else
				    logevent
					("Pageant's response not accepted");
			    } else {
				logevent
				    ("Pageant failed to answer challenge");
				sfree(ret);
			    }
			} else {
			    logevent("No reply received from Pageant");
			}
		    }
		    freebn(s->key.exponent);
		    freebn(s->key.modulus);
		    freebn(s->challenge);
		    if (s->authed)
			break;
		}
		sfree(s->response);
		if (s->publickey_blob && !s->tried_publickey)
		    logevent("Configured key file not in Pageant");
	    } else {
                logevent("Failed to get reply from Pageant");
            }
	    if (s->authed)
		break;
	}
	if (s->publickey_blob && !s->tried_publickey) {
	    /*
	     * Try public key authentication with the specified
	     * key file.
	     */
	    int got_passphrase; /* need not be kept over crReturn */
	    if (flags & FLAG_VERBOSE)
		c_write_str(ssh, "Trying public key authentication.\r\n");
	    logeventf(ssh, "Trying public key \"%s\"",
		      filename_to_str(&ssh->cfg.keyfile));
	    s->tried_publickey = 1;
	    got_passphrase = FALSE;
	    while (!got_passphrase) {
		/*
		 * Get a passphrase, if necessary.
		 */
		char *passphrase = NULL;    /* only written after crReturn */
		const char *error;
		if (!s->publickey_encrypted) {
		    if (flags & FLAG_VERBOSE)
			c_write_str(ssh, "No passphrase required.\r\n");
		    passphrase = NULL;
		} else {
		    int ret; /* need not be kept over crReturn */
		    s->cur_prompt = new_prompts(ssh->frontend);
		    s->cur_prompt->to_server = FALSE;
		    s->cur_prompt->name = dupstr("SSH key passphrase");
		    add_prompt(s->cur_prompt,
			       dupprintf("Passphrase for key \"%.100s\": ",
					 s->publickey_comment),
			       FALSE, SSH_MAX_PASSWORD_LEN);
		    ret = get_userpass_input(s->cur_prompt, NULL, 0);
		    while (ret < 0) {
			ssh->send_ok = 1;
			crWaitUntil(!pktin);
			ret = get_userpass_input(s->cur_prompt, in, inlen);
			ssh->send_ok = 0;
		    }
		    if (!ret) {
			/* Failed to get a passphrase. Terminate. */
			free_prompts(s->cur_prompt);
			ssh_disconnect(ssh, NULL, "Unable to authenticate",
				       0, TRUE);
			crStop(0);
		    }
		    passphrase = dupstr(s->cur_prompt->prompts[0]->result);
		    free_prompts(s->cur_prompt);
		}
		/*
		 * Try decrypting key with passphrase.
		 */
		ret = loadrsakey(&ssh->cfg.keyfile, &s->key, passphrase,
				 &error);
		if (passphrase) {
		    memset(passphrase, 0, strlen(passphrase));
		    sfree(passphrase);
		}
		if (ret == 1) {
		    /* Correct passphrase. */
		    got_passphrase = TRUE;
		} else if (ret == 0) {
		    c_write_str(ssh, "Couldn't load private key from ");
		    c_write_str(ssh, filename_to_str(&ssh->cfg.keyfile));
		    c_write_str(ssh, " (");
		    c_write_str(ssh, error);
		    c_write_str(ssh, ").\r\n");
		    got_passphrase = FALSE;
		    break;	       /* go and try something else */
		} else if (ret == -1) {
		    c_write_str(ssh, "Wrong passphrase.\r\n"); /* FIXME */
		    got_passphrase = FALSE;
		    /* and try again */
		} else {
		    assert(0 && "unexpected return from loadrsakey()");
		    got_passphrase = FALSE;   /* placate optimisers */
		}
	    }

	    if (got_passphrase) {

		/*
		 * Send a public key attempt.
		 */
		send_packet(ssh, SSH1_CMSG_AUTH_RSA,
			    PKT_BIGNUM, s->key.modulus, PKT_END);

		crWaitUntil(pktin);
		if (pktin->type == SSH1_SMSG_FAILURE) {
		    c_write_str(ssh, "Server refused our public key.\r\n");
		    continue;	       /* go and try something else */
		}
		if (pktin->type != SSH1_SMSG_AUTH_RSA_CHALLENGE) {
		    bombout(("Bizarre response to offer of public key"));
		    crStop(0);
		}

		{
		    int i;
		    unsigned char buffer[32];
		    Bignum challenge, response;

		    if ((challenge = ssh1_pkt_getmp(pktin)) == NULL) {
			bombout(("Server's RSA challenge was badly formatted"));
			crStop(0);
		    }
		    response = rsadecrypt(challenge, &s->key);
		    freebn(s->key.private_exponent);/* burn the evidence */

		    for (i = 0; i < 32; i++) {
			buffer[i] = bignum_byte(response, 31 - i);
		    }

		    MD5Init(&md5c);
		    MD5Update(&md5c, buffer, 32);
		    MD5Update(&md5c, s->session_id, 16);
		    MD5Final(buffer, &md5c);

		    send_packet(ssh, SSH1_CMSG_AUTH_RSA_RESPONSE,
				PKT_DATA, buffer, 16, PKT_END);

		    freebn(challenge);
		    freebn(response);
		}

		crWaitUntil(pktin);
		if (pktin->type == SSH1_SMSG_FAILURE) {
		    if (flags & FLAG_VERBOSE)
			c_write_str(ssh, "Failed to authenticate with"
				    " our public key.\r\n");
		    continue;	       /* go and try something else */
		} else if (pktin->type != SSH1_SMSG_SUCCESS) {
		    bombout(("Bizarre response to RSA authentication response"));
		    crStop(0);
		}

		break;		       /* we're through! */
	    }

	}

	/*
	 * Otherwise, try various forms of password-like authentication.
	 */
	s->cur_prompt = new_prompts(ssh->frontend);

	if (ssh->cfg.try_tis_auth &&
	    (s->supported_auths_mask & (1 << SSH1_AUTH_TIS)) &&
	    !s->tis_auth_refused) {
	    s->pwpkt_type = SSH1_CMSG_AUTH_TIS_RESPONSE;
	    logevent("Requested TIS authentication");
	    send_packet(ssh, SSH1_CMSG_AUTH_TIS, PKT_END);
	    crWaitUntil(pktin);
	    if (pktin->type != SSH1_SMSG_AUTH_TIS_CHALLENGE) {
		logevent("TIS authentication declined");
		if (flags & FLAG_INTERACTIVE)
		    c_write_str(ssh, "TIS authentication refused.\r\n");
		s->tis_auth_refused = 1;
		continue;
	    } else {
		char *challenge;
		int challengelen;
		char *instr_suf, *prompt;

		ssh_pkt_getstring(pktin, &challenge, &challengelen);
		if (!challenge) {
		    bombout(("TIS challenge packet was badly formed"));
		    crStop(0);
		}
		logevent("Received TIS challenge");
		s->cur_prompt->to_server = TRUE;
		s->cur_prompt->name = dupstr("SSH TIS authentication");
		/* Prompt heuristic comes from OpenSSH */
		if (memchr(challenge, '\n', challengelen)) {
		    instr_suf = dupstr("");
		    prompt = dupprintf("%.*s", challengelen, challenge);
		} else {
		    instr_suf = dupprintf("%.*s", challengelen, challenge);
		    prompt = dupstr("Response: ");
		}
		s->cur_prompt->instruction =
		    dupprintf("Using TIS authentication.%s%s",
			      (*instr_suf) ? "\n" : "",
			      instr_suf);
		s->cur_prompt->instr_reqd = TRUE;
		add_prompt(s->cur_prompt, prompt, FALSE, SSH_MAX_PASSWORD_LEN);
		sfree(instr_suf);
	    }
	}
	if (ssh->cfg.try_tis_auth &&
	    (s->supported_auths_mask & (1 << SSH1_AUTH_CCARD)) &&
	    !s->ccard_auth_refused) {
	    s->pwpkt_type = SSH1_CMSG_AUTH_CCARD_RESPONSE;
	    logevent("Requested CryptoCard authentication");
	    send_packet(ssh, SSH1_CMSG_AUTH_CCARD, PKT_END);
	    crWaitUntil(pktin);
	    if (pktin->type != SSH1_SMSG_AUTH_CCARD_CHALLENGE) {
		logevent("CryptoCard authentication declined");
		c_write_str(ssh, "CryptoCard authentication refused.\r\n");
		s->ccard_auth_refused = 1;
		continue;
	    } else {
		char *challenge;
		int challengelen;
		char *instr_suf, *prompt;

		ssh_pkt_getstring(pktin, &challenge, &challengelen);
		if (!challenge) {
		    bombout(("CryptoCard challenge packet was badly formed"));
		    crStop(0);
		}
		logevent("Received CryptoCard challenge");
		s->cur_prompt->to_server = TRUE;
		s->cur_prompt->name = dupstr("SSH CryptoCard authentication");
		s->cur_prompt->name_reqd = FALSE;
		/* Prompt heuristic comes from OpenSSH */
		if (memchr(challenge, '\n', challengelen)) {
		    instr_suf = dupstr("");
		    prompt = dupprintf("%.*s", challengelen, challenge);
		} else {
		    instr_suf = dupprintf("%.*s", challengelen, challenge);
		    prompt = dupstr("Response: ");
		}
		s->cur_prompt->instruction =
		    dupprintf("Using CryptoCard authentication.%s%s",
			      (*instr_suf) ? "\n" : "",
			      instr_suf);
		s->cur_prompt->instr_reqd = TRUE;
		add_prompt(s->cur_prompt, prompt, FALSE, SSH_MAX_PASSWORD_LEN);
		sfree(instr_suf);
	    }
	}
	if (s->pwpkt_type == SSH1_CMSG_AUTH_PASSWORD) {
	    if ((s->supported_auths_mask & (1 << SSH1_AUTH_PASSWORD)) == 0) {
		bombout(("No supported authentication methods available"));
		crStop(0);
	    }
	    s->cur_prompt->to_server = TRUE;
	    s->cur_prompt->name = dupstr("SSH password");
	    add_prompt(s->cur_prompt, dupprintf("%.90s@%.90s's password: ",
						s->username, ssh->savedhost),
		       FALSE, SSH_MAX_PASSWORD_LEN);
	}

	/*
	 * Show password prompt, having first obtained it via a TIS
	 * or CryptoCard exchange if we're doing TIS or CryptoCard
	 * authentication.
	 */
	{
	    int ret; /* need not be kept over crReturn */
	    ret = get_userpass_input(s->cur_prompt, NULL, 0);
	    while (ret < 0) {
		ssh->send_ok = 1;
		crWaitUntil(!pktin);
		ret = get_userpass_input(s->cur_prompt, in, inlen);
		ssh->send_ok = 0;
	    }
	    if (!ret) {
		/*
		 * Failed to get a password (for example
		 * because one was supplied on the command line
		 * which has already failed to work). Terminate.
		 */
		free_prompts(s->cur_prompt);
		ssh_disconnect(ssh, NULL, "Unable to authenticate", 0, TRUE);
		crStop(0);
	    }
	}

	if (s->pwpkt_type == SSH1_CMSG_AUTH_PASSWORD) {
	    /*
	     * Defence against traffic analysis: we send a
	     * whole bunch of packets containing strings of
	     * different lengths. One of these strings is the
	     * password, in a SSH1_CMSG_AUTH_PASSWORD packet.
	     * The others are all random data in
	     * SSH1_MSG_IGNORE packets. This way a passive
	     * listener can't tell which is the password, and
	     * hence can't deduce the password length.
	     * 
	     * Anybody with a password length greater than 16
	     * bytes is going to have enough entropy in their
	     * password that a listener won't find it _that_
	     * much help to know how long it is. So what we'll
	     * do is:
	     * 
	     *  - if password length < 16, we send 15 packets
	     *    containing string lengths 1 through 15
	     * 
	     *  - otherwise, we let N be the nearest multiple
	     *    of 8 below the password length, and send 8
	     *    packets containing string lengths N through
	     *    N+7. This won't obscure the order of
	     *    magnitude of the password length, but it will
	     *    introduce a bit of extra uncertainty.
	     * 
	     * A few servers can't deal with SSH1_MSG_IGNORE, at
	     * least in this context. For these servers, we need
	     * an alternative defence. We make use of the fact
	     * that the password is interpreted as a C string:
	     * so we can append a NUL, then some random data.
	     * 
	     * A few servers can deal with neither SSH1_MSG_IGNORE
	     * here _nor_ a padded password string.
	     * For these servers we are left with no defences
	     * against password length sniffing.
	     */
	    if (!(ssh->remote_bugs & BUG_CHOKES_ON_SSH1_IGNORE) &&
	        !(ssh->remote_bugs & BUG_NEEDS_SSH1_PLAIN_PASSWORD)) {
		/*
		 * The server can deal with SSH1_MSG_IGNORE, so
		 * we can use the primary defence.
		 */
		int bottom, top, pwlen, i;
		char *randomstr;

		pwlen = strlen(s->cur_prompt->prompts[0]->result);
		if (pwlen < 16) {
		    bottom = 0;    /* zero length passwords are OK! :-) */
		    top = 15;
		} else {
		    bottom = pwlen & ~7;
		    top = bottom + 7;
		}

		assert(pwlen >= bottom && pwlen <= top);

		randomstr = snewn(top + 1, char);

		for (i = bottom; i <= top; i++) {
		    if (i == pwlen) {
			defer_packet(ssh, s->pwpkt_type,
				     PKTT_PASSWORD, PKT_STR,
				     s->cur_prompt->prompts[0]->result,
				     PKTT_OTHER, PKT_END);
		    } else {
			for (j = 0; j < i; j++) {
			    do {
				randomstr[j] = random_byte();
			    } while (randomstr[j] == '\0');
			}
			randomstr[i] = '\0';
			defer_packet(ssh, SSH1_MSG_IGNORE,
				     PKT_STR, randomstr, PKT_END);
		    }
		}
		logevent("Sending password with camouflage packets");
		ssh_pkt_defersend(ssh);
		sfree(randomstr);
	    } 
	    else if (!(ssh->remote_bugs & BUG_NEEDS_SSH1_PLAIN_PASSWORD)) {
		/*
		 * The server can't deal with SSH1_MSG_IGNORE
		 * but can deal with padded passwords, so we
		 * can use the secondary defence.
		 */
		char string[64];
		char *ss;
		int len;

		len = strlen(s->cur_prompt->prompts[0]->result);
		if (len < sizeof(string)) {
		    ss = string;
		    strcpy(string, s->cur_prompt->prompts[0]->result);
		    len++;	       /* cover the zero byte */
		    while (len < sizeof(string)) {
			string[len++] = (char) random_byte();
		    }
		} else {
		    ss = s->cur_prompt->prompts[0]->result;
		}
		logevent("Sending length-padded password");
		send_packet(ssh, s->pwpkt_type, PKTT_PASSWORD,
			    PKT_INT, len, PKT_DATA, ss, len,
			    PKTT_OTHER, PKT_END);
	    } else {
		/*
		 * The server is believed unable to cope with
		 * any of our password camouflage methods.
		 */
		int len;
		len = strlen(s->cur_prompt->prompts[0]->result);
		logevent("Sending unpadded password");
		send_packet(ssh, s->pwpkt_type,
			    PKTT_PASSWORD, PKT_INT, len,
			    PKT_DATA, s->cur_prompt->prompts[0]->result, len,
			    PKTT_OTHER, PKT_END);
	    }
	} else {
	    send_packet(ssh, s->pwpkt_type, PKTT_PASSWORD,
			PKT_STR, s->cur_prompt->prompts[0]->result,
			PKTT_OTHER, PKT_END);
	}
	logevent("Sent password");
	free_prompts(s->cur_prompt);
	crWaitUntil(pktin);
	if (pktin->type == SSH1_SMSG_FAILURE) {
	    if (flags & FLAG_VERBOSE)
		c_write_str(ssh, "Access denied\r\n");
	    logevent("Authentication refused");
	} else if (pktin->type != SSH1_SMSG_SUCCESS) {
	    bombout(("Strange packet received, type %d", pktin->type));
	    crStop(0);
	}
    }

    /* Clear up */
    if (s->publickey_blob) {
	sfree(s->publickey_blob);
	sfree(s->publickey_comment);
    }

    logevent("Authentication successful");

    crFinish(1);
}

void sshfwd_close(struct ssh_channel *c)
{
    Ssh ssh = c->ssh;

    if (ssh->state == SSH_STATE_CLOSED)
	return;

    if (!c->closes) {
	/*
	 * If halfopen is true, we have sent
	 * CHANNEL_OPEN for this channel, but it hasn't even been
	 * acknowledged by the server. So we must set a close flag
	 * on it now, and then when the server acks the channel
	 * open, we can close it then.
	 */
	if (!c->halfopen) {
	    if (ssh->version == 1) {
		send_packet(ssh, SSH1_MSG_CHANNEL_CLOSE, PKT_INT, c->remoteid,
			    PKT_END);
		c->closes = 1;		       /* sent MSG_CLOSE */
	    } else {
		int bytes_to_send = bufchain_size(&c->v.v2.outbuffer);
		if (bytes_to_send > 0) {
		    /*
		     * If we still have unsent data in our outgoing
		     * buffer for this channel, we can't actually
		     * initiate a close operation yet or that data
		     * will be lost. Instead, set the pending_close
		     * flag so that when we do clear the buffer
		     * we'll start closing the channel.
		     */
		    char logmsg[160] = {'\0'};
		    sprintf(
			    logmsg,
			    "Forwarded port pending to be closed : "
			    "%d bytes remaining",
			    bytes_to_send);
		    logevent(logmsg);

		    c->pending_close = TRUE;
		} else {
		    /*
		     * No locally buffered data, so we can send the
		     * close message immediately.
		     */
		    struct Packet *pktout;
		    pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_CLOSE);
		    ssh2_pkt_adduint32(pktout, c->remoteid);
		    ssh2_pkt_send(ssh, pktout);
		    c->closes = 1;		       /* sent MSG_CLOSE */
		    logevent("Nothing left to send, closing channel");
		}
	    }
	}

	if (c->type == CHAN_X11) {
	    c->u.x11.s = NULL;
	    logevent("Forwarded X11 connection terminated");
	} else if (c->type == CHAN_SOCKDATA ||
		   c->type == CHAN_SOCKDATA_DORMANT) {
	    c->u.pfd.s = NULL;
	    logevent("Forwarded port closed");
	}
    }
}

int sshfwd_write(struct ssh_channel *c, char *buf, int len)
{
    Ssh ssh = c->ssh;

    if (ssh->state == SSH_STATE_CLOSED)
	return 0;

    if (ssh->version == 1) {
	send_packet(ssh, SSH1_MSG_CHANNEL_DATA,
		    PKT_INT, c->remoteid,
		    PKT_INT, len, PKTT_DATA, PKT_DATA, buf, len,
		    PKTT_OTHER, PKT_END);
	/*
	 * In SSH-1 we can return 0 here - implying that forwarded
	 * connections are never individually throttled - because
	 * the only circumstance that can cause throttling will be
	 * the whole SSH connection backing up, in which case
	 * _everything_ will be throttled as a whole.
	 */
	return 0;
    } else {
	ssh2_add_channel_data(c, buf, len);
	return ssh2_try_send(c);
    }
}

void sshfwd_unthrottle(struct ssh_channel *c, int bufsize)
{
    Ssh ssh = c->ssh;
    int buflimit;

    if (ssh->state == SSH_STATE_CLOSED)
	return;

    if (ssh->version == 1) {
	buflimit = SSH1_BUFFER_LIMIT;
    } else {
	buflimit = c->v.v2.locmaxwin;
	ssh2_set_window(c, bufsize < buflimit ? buflimit - bufsize : 0);
    }
    if (c->throttling_conn && bufsize <= buflimit) {
	c->throttling_conn = 0;
	ssh_throttle_conn(ssh, -1);
    }
}

static void ssh_queueing_handler(Ssh ssh, struct Packet *pktin)
{
    struct queued_handler *qh = ssh->qhead;

    assert(qh != NULL);

    assert(pktin->type == qh->msg1 || pktin->type == qh->msg2);

    if (qh->msg1 > 0) {
	assert(ssh->packet_dispatch[qh->msg1] == ssh_queueing_handler);
	ssh->packet_dispatch[qh->msg1] = NULL;
    }
    if (qh->msg2 > 0) {
	assert(ssh->packet_dispatch[qh->msg2] == ssh_queueing_handler);
	ssh->packet_dispatch[qh->msg2] = NULL;
    }

    if (qh->next) {
	ssh->qhead = qh->next;

	if (ssh->qhead->msg1 > 0) {
	    assert(ssh->packet_dispatch[ssh->qhead->msg1] == NULL);
	    ssh->packet_dispatch[ssh->qhead->msg1] = ssh_queueing_handler;
	}
	if (ssh->qhead->msg2 > 0) {
	    assert(ssh->packet_dispatch[ssh->qhead->msg2] == NULL);
	    ssh->packet_dispatch[ssh->qhead->msg2] = ssh_queueing_handler;
	}
    } else {
	ssh->qhead = ssh->qtail = NULL;
	ssh->packet_dispatch[pktin->type] = NULL;
    }

    qh->handler(ssh, pktin, qh->ctx);

    sfree(qh);
}

static void ssh_queue_handler(Ssh ssh, int msg1, int msg2,
			      chandler_fn_t handler, void *ctx)
{
    struct queued_handler *qh;

    qh = snew(struct queued_handler);
    qh->msg1 = msg1;
    qh->msg2 = msg2;
    qh->handler = handler;
    qh->ctx = ctx;
    qh->next = NULL;

    if (ssh->qtail == NULL) {
	ssh->qhead = qh;

	if (qh->msg1 > 0) {
	    assert(ssh->packet_dispatch[qh->msg1] == NULL);
	    ssh->packet_dispatch[qh->msg1] = ssh_queueing_handler;
	}
	if (qh->msg2 > 0) {
	    assert(ssh->packet_dispatch[qh->msg2] == NULL);
	    ssh->packet_dispatch[qh->msg2] = ssh_queueing_handler;
	}
    } else {
	ssh->qtail->next = qh;
    }
    ssh->qtail = qh;
}

static void ssh_rportfwd_succfail(Ssh ssh, struct Packet *pktin, void *ctx)
{
    struct ssh_rportfwd *rpf, *pf = (struct ssh_rportfwd *)ctx;

    if (pktin->type == (ssh->version == 1 ? SSH1_SMSG_SUCCESS :
			SSH2_MSG_REQUEST_SUCCESS)) {
	logeventf(ssh, "Remote port forwarding from %s enabled",
		  pf->sportdesc);
    } else {
	logeventf(ssh, "Remote port forwarding from %s refused",
		  pf->sportdesc);

	rpf = del234(ssh->rportfwds, pf);
	assert(rpf == pf);
	pf->pfrec->remote = NULL;
	free_rportfwd(pf);
    }
}

static void ssh_setup_portfwd(Ssh ssh, const Config *cfg)
{
    const char *portfwd_strptr = cfg->portfwd;
    struct ssh_portfwd *epf;
    int i;

    if (!ssh->portfwds) {
	ssh->portfwds = newtree234(ssh_portcmp);
    } else {
	/*
	 * Go through the existing port forwardings and tag them
	 * with status==DESTROY. Any that we want to keep will be
	 * re-enabled (status==KEEP) as we go through the
	 * configuration and find out which bits are the same as
	 * they were before.
	 */
	struct ssh_portfwd *epf;
	int i;
	for (i = 0; (epf = index234(ssh->portfwds, i)) != NULL; i++)
	    epf->status = DESTROY;
    }

    while (*portfwd_strptr) {
	char address_family, type;
	int sport,dport,sserv,dserv;
	char sports[256], dports[256], saddr[256], host[256];
	int n;

	address_family = 'A';
	type = 'L';
	if (*portfwd_strptr == 'A' ||
	    *portfwd_strptr == '4' ||
	    *portfwd_strptr == '6')
	    address_family = *portfwd_strptr++;
	if (*portfwd_strptr == 'L' ||
	    *portfwd_strptr == 'R' ||
	    *portfwd_strptr == 'D')
	    type = *portfwd_strptr++;

	saddr[0] = '\0';

	n = 0;
	while (*portfwd_strptr && *portfwd_strptr != '\t') {
	    if (*portfwd_strptr == ':') {
		/*
		 * We've seen a colon in the middle of the
		 * source port number. This means that
		 * everything we've seen until now is the
		 * source _address_, so we'll move it into
		 * saddr and start sports from the beginning
		 * again.
		 */
		portfwd_strptr++;
		sports[n] = '\0';
		if (ssh->version == 1 && type == 'R') {
		    logeventf(ssh, "SSH-1 cannot handle remote source address "
			      "spec \"%s\"; ignoring", sports);
		} else
		    strcpy(saddr, sports);
		n = 0;
	    }
	    if (n < lenof(sports)-1) sports[n++] = *portfwd_strptr++;
	}
	sports[n] = 0;
	if (type != 'D') {
	    if (*portfwd_strptr == '\t')
		portfwd_strptr++;
	    n = 0;
	    while (*portfwd_strptr && *portfwd_strptr != ':') {
		if (n < lenof(host)-1) host[n++] = *portfwd_strptr++;
	    }
	    host[n] = 0;
	    if (*portfwd_strptr == ':')
		portfwd_strptr++;
	    n = 0;
	    while (*portfwd_strptr) {
		if (n < lenof(dports)-1) dports[n++] = *portfwd_strptr++;
	    }
	    dports[n] = 0;
	    portfwd_strptr++;
	    dport = atoi(dports);
	    dserv = 0;
	    if (dport == 0) {
		dserv = 1;
		dport = net_service_lookup(dports);
		if (!dport) {
		    logeventf(ssh, "Service lookup failed for destination"
			      " port \"%s\"", dports);
		}
	    }
	} else {
	    while (*portfwd_strptr) portfwd_strptr++;
	    host[0] = 0;
	    dports[0] = 0;
	    dport = dserv = -1;
	    portfwd_strptr++;	       /* eat the NUL and move to next one */
	}
	sport = atoi(sports);
	sserv = 0;
	if (sport == 0) {
	    sserv = 1;
	    sport = net_service_lookup(sports);
	    if (!sport) {
		logeventf(ssh, "Service lookup failed for source"
			  " port \"%s\"", sports);
	    }
	}
	if (sport && dport) {
	    /* Set up a description of the source port. */
	    struct ssh_portfwd *pfrec, *epfrec;

	    pfrec = snew(struct ssh_portfwd);
	    pfrec->type = type;
	    pfrec->saddr = *saddr ? dupstr(saddr) : NULL;
	    pfrec->sserv = sserv ? dupstr(sports) : NULL;
	    pfrec->sport = sport;
	    pfrec->daddr = *host ? dupstr(host) : NULL;
	    pfrec->dserv = dserv ? dupstr(dports) : NULL;
	    pfrec->dport = dport;
	    pfrec->local = NULL;
	    pfrec->remote = NULL;
	    pfrec->addressfamily = (address_family == '4' ? ADDRTYPE_IPV4 :
				    address_family == '6' ? ADDRTYPE_IPV6 :
				    ADDRTYPE_UNSPEC);

	    epfrec = add234(ssh->portfwds, pfrec);
	    if (epfrec != pfrec) {
		if (epfrec->status == DESTROY) {
		    /*
		     * We already have a port forwarding up and running
		     * with precisely these parameters. Hence, no need
		     * to do anything; simply re-tag the existing one
		     * as KEEP.
		     */
		    epfrec->status = KEEP;
		}
		/*
		 * Anything else indicates that there was a duplicate
		 * in our input, which we'll silently ignore.
		 */
		free_portfwd(pfrec);
	    } else {
		pfrec->status = CREATE;
	    }
	}
    }

    /*
     * Now go through and destroy any port forwardings which were
     * not re-enabled.
     */
    for (i = 0; (epf = index234(ssh->portfwds, i)) != NULL; i++)
	if (epf->status == DESTROY) {
	    char *message;

	    message = dupprintf("%s port forwarding from %s%s%d",
				epf->type == 'L' ? "local" :
				epf->type == 'R' ? "remote" : "dynamic",
				epf->saddr ? epf->saddr : "",
				epf->saddr ? ":" : "",
				epf->sport);

	    if (epf->type != 'D') {
		char *msg2 = dupprintf("%s to %s:%d", message,
				       epf->daddr, epf->dport);
		sfree(message);
		message = msg2;
	    }

	    logeventf(ssh, "Cancelling %s", message);
	    sfree(message);

	    /* epf->remote or epf->local may be NULL if setting up a
	     * forwarding failed. */
	    if (epf->remote) {
		struct ssh_rportfwd *rpf = epf->remote;
		struct Packet *pktout;

		/*
		 * Cancel the port forwarding at the server
		 * end.
		 */
		if (ssh->version == 1) {
		    /*
		     * We cannot cancel listening ports on the
		     * server side in SSH-1! There's no message
		     * to support it. Instead, we simply remove
		     * the rportfwd record from the local end
		     * so that any connections the server tries
		     * to make on it are rejected.
		     */
		} else {
		    pktout = ssh2_pkt_init(SSH2_MSG_GLOBAL_REQUEST);
		    ssh2_pkt_addstring(pktout, "cancel-tcpip-forward");
		    ssh2_pkt_addbool(pktout, 0);/* _don't_ want reply */
		    if (epf->saddr) {
			ssh2_pkt_addstring(pktout, epf->saddr);
		    } else if (ssh->cfg.rport_acceptall) {
			/* XXX: ssh->cfg.rport_acceptall may not represent
			 * what was used to open the original connection,
			 * since it's reconfigurable. */
			ssh2_pkt_addstring(pktout, "0.0.0.0");
		    } else {
			ssh2_pkt_addstring(pktout, "127.0.0.1");
		    }
		    ssh2_pkt_adduint32(pktout, epf->sport);
		    ssh2_pkt_send(ssh, pktout);
		}

		del234(ssh->rportfwds, rpf);
		free_rportfwd(rpf);
	    } else if (epf->local) {
		pfd_terminate(epf->local);
	    }

	    delpos234(ssh->portfwds, i);
	    free_portfwd(epf);
	    i--;		       /* so we don't skip one in the list */
	}

    /*
     * And finally, set up any new port forwardings (status==CREATE).
     */
    for (i = 0; (epf = index234(ssh->portfwds, i)) != NULL; i++)
	if (epf->status == CREATE) {
	    char *sportdesc, *dportdesc;
	    sportdesc = dupprintf("%s%s%s%s%d%s",
				  epf->saddr ? epf->saddr : "",
				  epf->saddr ? ":" : "",
				  epf->sserv ? epf->sserv : "",
				  epf->sserv ? "(" : "",
				  epf->sport,
				  epf->sserv ? ")" : "");
	    if (epf->type == 'D') {
		dportdesc = NULL;
	    } else {
		dportdesc = dupprintf("%s:%s%s%d%s",
				      epf->daddr,
				      epf->dserv ? epf->dserv : "",
				      epf->dserv ? "(" : "",
				      epf->dport,
				      epf->dserv ? ")" : "");
	    }

	    if (epf->type == 'L') {
		const char *err = pfd_addforward(epf->daddr, epf->dport,
						 epf->saddr, epf->sport,
						 ssh, cfg,
						 &epf->local,
						 epf->addressfamily);

		logeventf(ssh, "Local %sport %s forwarding to %s%s%s",
			  epf->addressfamily == ADDRTYPE_IPV4 ? "IPv4 " :
			  epf->addressfamily == ADDRTYPE_IPV6 ? "IPv6 " : "",
			  sportdesc, dportdesc,
			  err ? " failed: " : "", err ? err : "");
	    } else if (epf->type == 'D') {
		const char *err = pfd_addforward(NULL, -1,
						 epf->saddr, epf->sport,
						 ssh, cfg,
						 &epf->local,
						 epf->addressfamily);

		logeventf(ssh, "Local %sport %s SOCKS dynamic forwarding%s%s",
			  epf->addressfamily == ADDRTYPE_IPV4 ? "IPv4 " :
			  epf->addressfamily == ADDRTYPE_IPV6 ? "IPv6 " : "",
			  sportdesc,
			  err ? " failed: " : "", err ? err : "");
	    } else {
		struct ssh_rportfwd *pf;

		/*
		 * Ensure the remote port forwardings tree exists.
		 */
		if (!ssh->rportfwds) {
		    if (ssh->version == 1)
			ssh->rportfwds = newtree234(ssh_rportcmp_ssh1);
		    else
			ssh->rportfwds = newtree234(ssh_rportcmp_ssh2);
		}

		pf = snew(struct ssh_rportfwd);
		strncpy(pf->dhost, epf->daddr, lenof(pf->dhost)-1);
		pf->dhost[lenof(pf->dhost)-1] = '\0';
		pf->dport = epf->dport;
		pf->sport = epf->sport;
		if (add234(ssh->rportfwds, pf) != pf) {
		    logeventf(ssh, "Duplicate remote port forwarding to %s:%d",
			      epf->daddr, epf->dport);
		    sfree(pf);
		} else {
		    logeventf(ssh, "Requesting remote port %s"
			      " forward to %s", sportdesc, dportdesc);

		    pf->sportdesc = sportdesc;
		    sportdesc = NULL;
		    epf->remote = pf;
		    pf->pfrec = epf;

		    if (ssh->version == 1) {
			send_packet(ssh, SSH1_CMSG_PORT_FORWARD_REQUEST,
				    PKT_INT, epf->sport,
				    PKT_STR, epf->daddr,
				    PKT_INT, epf->dport,
				    PKT_END);
			ssh_queue_handler(ssh, SSH1_SMSG_SUCCESS,
					  SSH1_SMSG_FAILURE,
					  ssh_rportfwd_succfail, pf);
		    } else {
			struct Packet *pktout;
			pktout = ssh2_pkt_init(SSH2_MSG_GLOBAL_REQUEST);
			ssh2_pkt_addstring(pktout, "tcpip-forward");
			ssh2_pkt_addbool(pktout, 1);/* want reply */
			if (epf->saddr) {
			    ssh2_pkt_addstring(pktout, epf->saddr);
			} else if (cfg->rport_acceptall) {
			    ssh2_pkt_addstring(pktout, "0.0.0.0");
			} else {
			    ssh2_pkt_addstring(pktout, "127.0.0.1");
			}
			ssh2_pkt_adduint32(pktout, epf->sport);
			ssh2_pkt_send(ssh, pktout);

			ssh_queue_handler(ssh, SSH2_MSG_REQUEST_SUCCESS,
					  SSH2_MSG_REQUEST_FAILURE,
					  ssh_rportfwd_succfail, pf);
		    }
		}
	    }
	    sfree(sportdesc);
	    sfree(dportdesc);
	}
}

static void ssh1_smsg_stdout_stderr_data(Ssh ssh, struct Packet *pktin)
{
    char *string;
    int stringlen, bufsize;

    ssh_pkt_getstring(pktin, &string, &stringlen);
    if (string == NULL) {
	bombout(("Incoming terminal data packet was badly formed"));
	return;
    }

    bufsize = from_backend(ssh->frontend, pktin->type == SSH1_SMSG_STDERR_DATA,
			   string, stringlen);
    if (!ssh->v1_stdout_throttling && bufsize > SSH1_BUFFER_LIMIT) {
	ssh->v1_stdout_throttling = 1;
	ssh_throttle_conn(ssh, +1);
    }
}

static void ssh1_smsg_x11_open(Ssh ssh, struct Packet *pktin)
{
    /* Remote side is trying to open a channel to talk to our
     * X-Server. Give them back a local channel number. */
    struct ssh_channel *c;
    int remoteid = ssh_pkt_getuint32(pktin);

    logevent("Received X11 connect request");
    /* Refuse if X11 forwarding is disabled. */
    if (!ssh->X11_fwd_enabled) {
	send_packet(ssh, SSH1_MSG_CHANNEL_OPEN_FAILURE,
		    PKT_INT, remoteid, PKT_END);
	logevent("Rejected X11 connect request");
    } else {
	c = snew(struct ssh_channel);
	c->ssh = ssh;

	if (x11_init(&c->u.x11.s, ssh->x11disp, c,
		     NULL, -1, &ssh->cfg) != NULL) {
	    logevent("Opening X11 forward connection failed");
	    sfree(c);
	    send_packet(ssh, SSH1_MSG_CHANNEL_OPEN_FAILURE,
			PKT_INT, remoteid, PKT_END);
	} else {
	    logevent
		("Opening X11 forward connection succeeded");
	    c->remoteid = remoteid;
	    c->halfopen = FALSE;
	    c->localid = alloc_channel_id(ssh);
	    c->closes = 0;
	    c->pending_close = FALSE;
	    c->throttling_conn = 0;
	    c->type = CHAN_X11;	/* identify channel type */
	    add234(ssh->channels, c);
	    send_packet(ssh, SSH1_MSG_CHANNEL_OPEN_CONFIRMATION,
			PKT_INT, c->remoteid, PKT_INT,
			c->localid, PKT_END);
	    logevent("Opened X11 forward channel");
	}
    }
}

static void ssh1_smsg_agent_open(Ssh ssh, struct Packet *pktin)
{
    /* Remote side is trying to open a channel to talk to our
     * agent. Give them back a local channel number. */
    struct ssh_channel *c;
    int remoteid = ssh_pkt_getuint32(pktin);

    /* Refuse if agent forwarding is disabled. */
    if (!ssh->agentfwd_enabled) {
	send_packet(ssh, SSH1_MSG_CHANNEL_OPEN_FAILURE,
		    PKT_INT, remoteid, PKT_END);
    } else {
	c = snew(struct ssh_channel);
	c->ssh = ssh;
	c->remoteid = remoteid;
	c->halfopen = FALSE;
	c->localid = alloc_channel_id(ssh);
	c->closes = 0;
	c->pending_close = FALSE;
	c->throttling_conn = 0;
	c->type = CHAN_AGENT;	/* identify channel type */
	c->u.a.lensofar = 0;
	add234(ssh->channels, c);
	send_packet(ssh, SSH1_MSG_CHANNEL_OPEN_CONFIRMATION,
		    PKT_INT, c->remoteid, PKT_INT, c->localid,
		    PKT_END);
    }
}

static void ssh1_msg_port_open(Ssh ssh, struct Packet *pktin)
{
    /* Remote side is trying to open a channel to talk to a
     * forwarded port. Give them back a local channel number. */
    struct ssh_channel *c;
    struct ssh_rportfwd pf, *pfp;
    int remoteid;
    int hostsize, port;
    char *host;
    const char *e;
    c = snew(struct ssh_channel);
    c->ssh = ssh;

    remoteid = ssh_pkt_getuint32(pktin);
    ssh_pkt_getstring(pktin, &host, &hostsize);
    port = ssh_pkt_getuint32(pktin);

    if (hostsize >= lenof(pf.dhost))
	hostsize = lenof(pf.dhost)-1;
    memcpy(pf.dhost, host, hostsize);
    pf.dhost[hostsize] = '\0';
    pf.dport = port;
    pfp = find234(ssh->rportfwds, &pf, NULL);

    if (pfp == NULL) {
	logeventf(ssh, "Rejected remote port open request for %s:%d",
		  pf.dhost, port);
	send_packet(ssh, SSH1_MSG_CHANNEL_OPEN_FAILURE,
		    PKT_INT, remoteid, PKT_END);
    } else {
	logeventf(ssh, "Received remote port open request for %s:%d",
		  pf.dhost, port);
	e = pfd_newconnect(&c->u.pfd.s, pf.dhost, port,
			   c, &ssh->cfg, pfp->pfrec->addressfamily);
	if (e != NULL) {
	    logeventf(ssh, "Port open failed: %s", e);
	    sfree(c);
	    send_packet(ssh, SSH1_MSG_CHANNEL_OPEN_FAILURE,
			PKT_INT, remoteid, PKT_END);
	} else {
	    c->remoteid = remoteid;
	    c->halfopen = FALSE;
	    c->localid = alloc_channel_id(ssh);
	    c->closes = 0;
	    c->pending_close = FALSE;
	    c->throttling_conn = 0;
	    c->type = CHAN_SOCKDATA;	/* identify channel type */
	    add234(ssh->channels, c);
	    send_packet(ssh, SSH1_MSG_CHANNEL_OPEN_CONFIRMATION,
			PKT_INT, c->remoteid, PKT_INT,
			c->localid, PKT_END);
	    logevent("Forwarded port opened successfully");
	}
    }
}

static void ssh1_msg_channel_open_confirmation(Ssh ssh, struct Packet *pktin)
{
    unsigned int remoteid = ssh_pkt_getuint32(pktin);
    unsigned int localid = ssh_pkt_getuint32(pktin);
    struct ssh_channel *c;

    c = find234(ssh->channels, &remoteid, ssh_channelfind);
    if (c && c->type == CHAN_SOCKDATA_DORMANT) {
	c->remoteid = localid;
	c->halfopen = FALSE;
	c->type = CHAN_SOCKDATA;
	c->throttling_conn = 0;
	pfd_confirm(c->u.pfd.s);
    }

    if (c && c->closes) {
	/*
	 * We have a pending close on this channel,
	 * which we decided on before the server acked
	 * the channel open. So now we know the
	 * remoteid, we can close it again.
	 */
	send_packet(ssh, SSH1_MSG_CHANNEL_CLOSE,
		    PKT_INT, c->remoteid, PKT_END);
    }
}

static void ssh1_msg_channel_open_failure(Ssh ssh, struct Packet *pktin)
{
    unsigned int remoteid = ssh_pkt_getuint32(pktin);
    struct ssh_channel *c;

    c = find234(ssh->channels, &remoteid, ssh_channelfind);
    if (c && c->type == CHAN_SOCKDATA_DORMANT) {
	logevent("Forwarded connection refused by server");
	pfd_close(c->u.pfd.s);
	del234(ssh->channels, c);
	sfree(c);
    }
}

static void ssh1_msg_channel_close(Ssh ssh, struct Packet *pktin)
{
    /* Remote side closes a channel. */
    unsigned i = ssh_pkt_getuint32(pktin);
    struct ssh_channel *c;
    c = find234(ssh->channels, &i, ssh_channelfind);
    if (c && !c->halfopen) {
	int closetype;
	closetype =
	    (pktin->type == SSH1_MSG_CHANNEL_CLOSE ? 1 : 2);

	if ((c->closes == 0) && (c->type == CHAN_X11)) {
	    logevent("Forwarded X11 connection terminated");
	    assert(c->u.x11.s != NULL);
	    x11_close(c->u.x11.s);
	    c->u.x11.s = NULL;
	}
	if ((c->closes == 0) && (c->type == CHAN_SOCKDATA)) {
	    logevent("Forwarded port closed");
	    assert(c->u.pfd.s != NULL);
	    pfd_close(c->u.pfd.s);
	    c->u.pfd.s = NULL;
	}

	c->closes |= (closetype << 2);   /* seen this message */
	if (!(c->closes & closetype)) {
	    send_packet(ssh, pktin->type, PKT_INT, c->remoteid,
			PKT_END);
	    c->closes |= closetype;      /* sent it too */
	}

	if (c->closes == 15) {
	    del234(ssh->channels, c);
	    sfree(c);
	}
    } else {
	bombout(("Received CHANNEL_CLOSE%s for %s channel %d\n",
		 pktin->type == SSH1_MSG_CHANNEL_CLOSE ? "" :
		 "_CONFIRMATION", c ? "half-open" : "nonexistent",
		 i));
    }
}

static void ssh1_msg_channel_data(Ssh ssh, struct Packet *pktin)
{
    /* Data sent down one of our channels. */
    int i = ssh_pkt_getuint32(pktin);
    char *p;
    int len;
    struct ssh_channel *c;

    ssh_pkt_getstring(pktin, &p, &len);

    c = find234(ssh->channels, &i, ssh_channelfind);
    if (c) {
	int bufsize = 0;
	switch (c->type) {
	  case CHAN_X11:
	    bufsize = x11_send(c->u.x11.s, p, len);
	    break;
	  case CHAN_SOCKDATA:
	    bufsize = pfd_send(c->u.pfd.s, p, len);
	    break;
	  case CHAN_AGENT:
	    /* Data for an agent message. Buffer it. */
	    while (len > 0) {
		if (c->u.a.lensofar < 4) {
		    unsigned int l = min(4 - c->u.a.lensofar, (unsigned)len);
		    memcpy(c->u.a.msglen + c->u.a.lensofar, p,
			   l);
		    p += l;
		    len -= l;
		    c->u.a.lensofar += l;
		}
		if (c->u.a.lensofar == 4) {
		    c->u.a.totallen =
			4 + GET_32BIT(c->u.a.msglen);
		    c->u.a.message = snewn(c->u.a.totallen,
					   unsigned char);
		    memcpy(c->u.a.message, c->u.a.msglen, 4);
		}
		if (c->u.a.lensofar >= 4 && len > 0) {
		    unsigned int l =
			min(c->u.a.totallen - c->u.a.lensofar,
			    (unsigned)len);
		    memcpy(c->u.a.message + c->u.a.lensofar, p,
			   l);
		    p += l;
		    len -= l;
		    c->u.a.lensofar += l;
		}
		if (c->u.a.lensofar == c->u.a.totallen) {
		    void *reply;
		    int replylen;
		    if (agent_query(c->u.a.message,
				    c->u.a.totallen,
				    &reply, &replylen,
				    ssh_agentf_callback, c))
			ssh_agentf_callback(c, reply, replylen);
		    sfree(c->u.a.message);
		    c->u.a.lensofar = 0;
		}
	    }
	    bufsize = 0;   /* agent channels never back up */
	    break;
	}
	if (!c->throttling_conn && bufsize > SSH1_BUFFER_LIMIT) {
	    c->throttling_conn = 1;
	    ssh_throttle_conn(ssh, +1);
	}
    }
}

static void ssh1_smsg_exit_status(Ssh ssh, struct Packet *pktin)
{
    ssh->exitcode = ssh_pkt_getuint32(pktin);
    logeventf(ssh, "Server sent command exit status %d", ssh->exitcode);
    send_packet(ssh, SSH1_CMSG_EXIT_CONFIRMATION, PKT_END);
    /*
     * In case `helpful' firewalls or proxies tack
     * extra human-readable text on the end of the
     * session which we might mistake for another
     * encrypted packet, we close the session once
     * we've sent EXIT_CONFIRMATION.
     */
    ssh_disconnect(ssh, NULL, NULL, 0, TRUE);
}

/* Helper function to deal with sending tty modes for REQUEST_PTY */
static void ssh1_send_ttymode(void *data, char *mode, char *val)
{
    struct Packet *pktout = (struct Packet *)data;
    int i = 0;
    unsigned int arg = 0;
    while (strcmp(mode, ssh_ttymodes[i].mode) != 0) i++;
    if (i == lenof(ssh_ttymodes)) return;
    switch (ssh_ttymodes[i].type) {
      case TTY_OP_CHAR:
	arg = ssh_tty_parse_specchar(val);
	break;
      case TTY_OP_BOOL:
	arg = ssh_tty_parse_boolean(val);
	break;
    }
    ssh2_pkt_addbyte(pktout, ssh_ttymodes[i].opcode);
    ssh2_pkt_addbyte(pktout, arg);
}


static void do_ssh1_connection(Ssh ssh, unsigned char *in, int inlen,
			       struct Packet *pktin)
{
    crBegin(ssh->do_ssh1_connection_crstate);

    ssh->packet_dispatch[SSH1_SMSG_STDOUT_DATA] = 
	ssh->packet_dispatch[SSH1_SMSG_STDERR_DATA] =
	ssh1_smsg_stdout_stderr_data;

    ssh->packet_dispatch[SSH1_MSG_CHANNEL_OPEN_CONFIRMATION] =
	ssh1_msg_channel_open_confirmation;
    ssh->packet_dispatch[SSH1_MSG_CHANNEL_OPEN_FAILURE] =
	ssh1_msg_channel_open_failure;
    ssh->packet_dispatch[SSH1_MSG_CHANNEL_CLOSE] =
	ssh->packet_dispatch[SSH1_MSG_CHANNEL_CLOSE_CONFIRMATION] =
	ssh1_msg_channel_close;
    ssh->packet_dispatch[SSH1_MSG_CHANNEL_DATA] = ssh1_msg_channel_data;
    ssh->packet_dispatch[SSH1_SMSG_EXIT_STATUS] = ssh1_smsg_exit_status;

    if (ssh->cfg.agentfwd && agent_exists()) {
	logevent("Requesting agent forwarding");
	send_packet(ssh, SSH1_CMSG_AGENT_REQUEST_FORWARDING, PKT_END);
	do {
	    crReturnV;
	} while (!pktin);
	if (pktin->type != SSH1_SMSG_SUCCESS
	    && pktin->type != SSH1_SMSG_FAILURE) {
	    bombout(("Protocol confusion"));
	    crStopV;
	} else if (pktin->type == SSH1_SMSG_FAILURE) {
	    logevent("Agent forwarding refused");
	} else {
	    logevent("Agent forwarding enabled");
	    ssh->agentfwd_enabled = TRUE;
	    ssh->packet_dispatch[SSH1_SMSG_AGENT_OPEN] = ssh1_smsg_agent_open;
	}
    }

    if (ssh->cfg.x11_forward &&
	(ssh->x11disp = x11_setup_display(ssh->cfg.x11_display,
					  ssh->cfg.x11_auth, &ssh->cfg))) {
	logevent("Requesting X11 forwarding");
	/*
	 * Note that while we blank the X authentication data here, we don't
	 * take any special action to blank the start of an X11 channel,
	 * so using MIT-MAGIC-COOKIE-1 and actually opening an X connection
	 * without having session blanking enabled is likely to leak your
	 * cookie into the log.
	 */
	if (ssh->v1_local_protoflags & SSH1_PROTOFLAG_SCREEN_NUMBER) {
	    send_packet(ssh, SSH1_CMSG_X11_REQUEST_FORWARDING,
			PKT_STR, ssh->x11disp->remoteauthprotoname,
			PKTT_PASSWORD,
			PKT_STR, ssh->x11disp->remoteauthdatastring,
			PKTT_OTHER,
			PKT_INT, ssh->x11disp->screennum,
			PKT_END);
	} else {
	    send_packet(ssh, SSH1_CMSG_X11_REQUEST_FORWARDING,
			PKT_STR, ssh->x11disp->remoteauthprotoname,
			PKTT_PASSWORD,
			PKT_STR, ssh->x11disp->remoteauthdatastring,
			PKTT_OTHER,
			PKT_END);
	}
	do {
	    crReturnV;
	} while (!pktin);
	if (pktin->type != SSH1_SMSG_SUCCESS
	    && pktin->type != SSH1_SMSG_FAILURE) {
	    bombout(("Protocol confusion"));
	    crStopV;
	} else if (pktin->type == SSH1_SMSG_FAILURE) {
	    logevent("X11 forwarding refused");
	} else {
	    logevent("X11 forwarding enabled");
	    ssh->X11_fwd_enabled = TRUE;
	    ssh->packet_dispatch[SSH1_SMSG_X11_OPEN] = ssh1_smsg_x11_open;
	}
    }

    ssh_setup_portfwd(ssh, &ssh->cfg);
    ssh->packet_dispatch[SSH1_MSG_PORT_OPEN] = ssh1_msg_port_open;

    if (!ssh->cfg.nopty) {
	struct Packet *pkt;
	/* Unpick the terminal-speed string. */
	/* XXX perhaps we should allow no speeds to be sent. */
	ssh->ospeed = 38400; ssh->ispeed = 38400; /* last-resort defaults */
	sscanf(ssh->cfg.termspeed, "%d,%d", &ssh->ospeed, &ssh->ispeed);
	/* Send the pty request. */
	pkt = ssh1_pkt_init(SSH1_CMSG_REQUEST_PTY);
	ssh_pkt_addstring(pkt, ssh->cfg.termtype);
	ssh_pkt_adduint32(pkt, ssh->term_height);
	ssh_pkt_adduint32(pkt, ssh->term_width);
	ssh_pkt_adduint32(pkt, 0); /* width in pixels */
	ssh_pkt_adduint32(pkt, 0); /* height in pixels */
	parse_ttymodes(ssh, ssh->cfg.ttymodes,
		       ssh1_send_ttymode, (void *)pkt);
	ssh_pkt_addbyte(pkt, SSH1_TTY_OP_ISPEED);
	ssh_pkt_adduint32(pkt, ssh->ispeed);
	ssh_pkt_addbyte(pkt, SSH1_TTY_OP_OSPEED);
	ssh_pkt_adduint32(pkt, ssh->ospeed);
	ssh_pkt_addbyte(pkt, SSH_TTY_OP_END);
	s_wrpkt(ssh, pkt);
	ssh->state = SSH_STATE_INTERMED;
	do {
	    crReturnV;
	} while (!pktin);
	if (pktin->type != SSH1_SMSG_SUCCESS
	    && pktin->type != SSH1_SMSG_FAILURE) {
	    bombout(("Protocol confusion"));
	    crStopV;
	} else if (pktin->type == SSH1_SMSG_FAILURE) {
	    c_write_str(ssh, "Server refused to allocate pty\r\n");
	    ssh->editing = ssh->echoing = 1;
	}
	logeventf(ssh, "Allocated pty (ospeed %dbps, ispeed %dbps)",
		  ssh->ospeed, ssh->ispeed);
    } else {
	ssh->editing = ssh->echoing = 1;
    }

    if (ssh->cfg.compression) {
	send_packet(ssh, SSH1_CMSG_REQUEST_COMPRESSION, PKT_INT, 6, PKT_END);
	do {
	    crReturnV;
	} while (!pktin);
	if (pktin->type != SSH1_SMSG_SUCCESS
	    && pktin->type != SSH1_SMSG_FAILURE) {
	    bombout(("Protocol confusion"));
	    crStopV;
	} else if (pktin->type == SSH1_SMSG_FAILURE) {
	    c_write_str(ssh, "Server refused to compress\r\n");
	}
	logevent("Started compression");
	ssh->v1_compressing = TRUE;
	ssh->cs_comp_ctx = zlib_compress_init();
	logevent("Initialised zlib (RFC1950) compression");
	ssh->sc_comp_ctx = zlib_decompress_init();
	logevent("Initialised zlib (RFC1950) decompression");
    }

    /*
     * Start the shell or command.
     * 
     * Special case: if the first-choice command is an SSH-2
     * subsystem (hence not usable here) and the second choice
     * exists, we fall straight back to that.
     */
    {
	char *cmd = ssh->cfg.remote_cmd_ptr;

	if (!cmd) cmd = ssh->cfg.remote_cmd;
	
	if (ssh->cfg.ssh_subsys && ssh->cfg.remote_cmd_ptr2) {
	    cmd = ssh->cfg.remote_cmd_ptr2;
	    ssh->fallback_cmd = TRUE;
	}
	if (*cmd)
	    send_packet(ssh, SSH1_CMSG_EXEC_CMD, PKT_STR, cmd, PKT_END);
	else
	    send_packet(ssh, SSH1_CMSG_EXEC_SHELL, PKT_END);
	logevent("Started session");
    }

    ssh->state = SSH_STATE_SESSION;
    if (ssh->size_needed)
	ssh_size(ssh, ssh->term_width, ssh->term_height);
    if (ssh->eof_needed)
	ssh_special(ssh, TS_EOF);

    if (ssh->ldisc)
	ldisc_send(ssh->ldisc, NULL, 0, 0);/* cause ldisc to notice changes */
    ssh->send_ok = 1;
    ssh->channels = newtree234(ssh_channelcmp);
    while (1) {

	/*
	 * By this point, most incoming packets are already being
	 * handled by the dispatch table, and we need only pay
	 * attention to the unusual ones.
	 */

	crReturnV;
	if (pktin) {
	    if (pktin->type == SSH1_SMSG_SUCCESS) {
		/* may be from EXEC_SHELL on some servers */
	    } else if (pktin->type == SSH1_SMSG_FAILURE) {
		/* may be from EXEC_SHELL on some servers
		 * if no pty is available or in other odd cases. Ignore */
	    } else {
		bombout(("Strange packet received: type %d", pktin->type));
		crStopV;
	    }
	} else {
	    while (inlen > 0) {
		int len = min(inlen, 512);
		send_packet(ssh, SSH1_CMSG_STDIN_DATA,
			    PKT_INT, len,  PKTT_DATA, PKT_DATA, in, len,
			    PKTT_OTHER, PKT_END);
		in += len;
		inlen -= len;
	    }
	}
    }

    crFinishV;
}

/*
 * Handle the top-level SSH-2 protocol.
 */
static void ssh1_msg_debug(Ssh ssh, struct Packet *pktin)
{
    char *msg;
    int msglen;

    ssh_pkt_getstring(pktin, &msg, &msglen);
    logeventf(ssh, "Remote debug message: %.*s", msglen, msg);
}

static void ssh1_msg_disconnect(Ssh ssh, struct Packet *pktin)
{
    /* log reason code in disconnect message */
    char *msg;
    int msglen;

    ssh_pkt_getstring(pktin, &msg, &msglen);
    bombout(("Server sent disconnect message:\n\"%.*s\"", msglen, msg));
}

static void ssh_msg_ignore(Ssh ssh, struct Packet *pktin)
{
    /* Do nothing, because we're ignoring it! Duhh. */
}

static void ssh1_protocol_setup(Ssh ssh)
{
    int i;

    /*
     * Most messages are handled by the coroutines.
     */
    for (i = 0; i < 256; i++)
	ssh->packet_dispatch[i] = NULL;

    /*
     * These special message types we install handlers for.
     */
    ssh->packet_dispatch[SSH1_MSG_DISCONNECT] = ssh1_msg_disconnect;
    ssh->packet_dispatch[SSH1_MSG_IGNORE] = ssh_msg_ignore;
    ssh->packet_dispatch[SSH1_MSG_DEBUG] = ssh1_msg_debug;
}

static void ssh1_protocol(Ssh ssh, void *vin, int inlen,
			  struct Packet *pktin)
{
    unsigned char *in=(unsigned char*)vin;
    if (ssh->state == SSH_STATE_CLOSED)
	return;

    if (pktin && ssh->packet_dispatch[pktin->type]) {
	ssh->packet_dispatch[pktin->type](ssh, pktin);
	return;
    }

    if (!ssh->protocol_initial_phase_done) {
	if (do_ssh1_login(ssh, in, inlen, pktin))
	    ssh->protocol_initial_phase_done = TRUE;
	else
	    return;
    }

    do_ssh1_connection(ssh, in, inlen, pktin);
}

/*
 * Utility routine for decoding comma-separated strings in KEXINIT.
 */
static int in_commasep_string(char *needle, char *haystack, int haylen)
{
    int needlen;
    if (!needle || !haystack)	       /* protect against null pointers */
	return 0;
    needlen = strlen(needle);
    while (1) {
	/*
	 * Is it at the start of the string?
	 */
	if (haylen >= needlen &&       /* haystack is long enough */
	    !memcmp(needle, haystack, needlen) &&	/* initial match */
	    (haylen == needlen || haystack[needlen] == ',')
	    /* either , or EOS follows */
	    )
	    return 1;
	/*
	 * If not, search for the next comma and resume after that.
	 * If no comma found, terminate.
	 */
	while (haylen > 0 && *haystack != ',')
	    haylen--, haystack++;
	if (haylen == 0)
	    return 0;
	haylen--, haystack++;	       /* skip over comma itself */
    }
}

/*
 * Similar routine for checking whether we have the first string in a list.
 */
static int first_in_commasep_string(char *needle, char *haystack, int haylen)
{
    int needlen;
    if (!needle || !haystack)	       /* protect against null pointers */
	return 0;
    needlen = strlen(needle);
    /*
     * Is it at the start of the string?
     */
    if (haylen >= needlen &&       /* haystack is long enough */
	!memcmp(needle, haystack, needlen) &&	/* initial match */
	(haylen == needlen || haystack[needlen] == ',')
	/* either , or EOS follows */
	)
	return 1;
    return 0;
}


/*
 * SSH-2 key creation method.
 * (Currently assumes 2 lots of any hash are sufficient to generate
 * keys/IVs for any cipher/MAC. SSH2_MKKEY_ITERS documents this assumption.)
 */
#define SSH2_MKKEY_ITERS (2)
static void ssh2_mkkey(Ssh ssh, Bignum K, unsigned char *H, char chr,
		       unsigned char *keyspace)
{
    const struct ssh_hash *h = ssh->kex->hash;
    void *s;
    /* First hlen bytes. */
    s = h->init();
    if (!(ssh->remote_bugs & BUG_SSH2_DERIVEKEY))
	hash_mpint(h, s, K);
    h->bytes(s, H, h->hlen);
    h->bytes(s, &chr, 1);
    h->bytes(s, ssh->v2_session_id, ssh->v2_session_id_len);
    h->final(s, keyspace);
    /* Next hlen bytes. */
    s = h->init();
    if (!(ssh->remote_bugs & BUG_SSH2_DERIVEKEY))
	hash_mpint(h, s, K);
    h->bytes(s, H, h->hlen);
    h->bytes(s, keyspace, h->hlen);
    h->final(s, keyspace + h->hlen);
}

/*
 * Handle the SSH-2 transport layer.
 */
static int do_ssh2_transport(Ssh ssh, void *vin, int inlen,
			     struct Packet *pktin)
{
    unsigned char *in = (unsigned char *)vin;
    struct do_ssh2_transport_state {
	int nbits, pbits, warn_kex, warn_cscipher, warn_sccipher;
	Bignum p, g, e, f, K;
	void *our_kexinit;
	int our_kexinitlen;
	int kex_init_value, kex_reply_value;
	const struct ssh_mac **maclist;
	int nmacs;
	const struct ssh2_cipher *cscipher_tobe;
	const struct ssh2_cipher *sccipher_tobe;
	const struct ssh_mac *csmac_tobe;
	const struct ssh_mac *scmac_tobe;
	const struct ssh_compress *cscomp_tobe;
	const struct ssh_compress *sccomp_tobe;
	char *hostkeydata, *sigdata, *rsakeydata, *keystr, *fingerprint;
	int hostkeylen, siglen, rsakeylen;
	void *hkey;		       /* actual host key */
	void *rsakey;		       /* for RSA kex */
	unsigned char exchange_hash[SSH2_KEX_MAX_HASH_LEN];
	int n_preferred_kex;
	const struct ssh_kexes *preferred_kex[KEX_MAX];
	int n_preferred_ciphers;
	const struct ssh2_ciphers *preferred_ciphers[CIPHER_MAX];
	const struct ssh_compress *preferred_comp;
	int userauth_succeeded;	    /* for delayed compression */
	int pending_compression;
	int got_session_id, activated_authconn;
	struct Packet *pktout;
        int dlgret;
	int guessok;
	int ignorepkt;
    };
    crState(do_ssh2_transport_state);

    crBegin(ssh->do_ssh2_transport_crstate);

    s->cscipher_tobe = s->sccipher_tobe = NULL;
    s->csmac_tobe = s->scmac_tobe = NULL;
    s->cscomp_tobe = s->sccomp_tobe = NULL;

    s->got_session_id = s->activated_authconn = FALSE;
    s->userauth_succeeded = FALSE;
    s->pending_compression = FALSE;

    /*
     * Be prepared to work around the buggy MAC problem.
     */
    if (ssh->remote_bugs & BUG_SSH2_HMAC)
	s->maclist = buggymacs, s->nmacs = lenof(buggymacs);
    else
	s->maclist = macs, s->nmacs = lenof(macs);

  begin_key_exchange:
    ssh->pkt_kctx = SSH2_PKTCTX_NOKEX;
    {
	int i, j, commalist_started;

	/*
	 * Set up the preferred key exchange. (NULL => warn below here)
	 */
	s->n_preferred_kex = 0;
	for (i = 0; i < KEX_MAX; i++) {
	    switch (ssh->cfg.ssh_kexlist[i]) {
	      case KEX_DHGEX:
		s->preferred_kex[s->n_preferred_kex++] =
		    &ssh_diffiehellman_gex;
		break;
	      case KEX_DHGROUP14:
		s->preferred_kex[s->n_preferred_kex++] =
		    &ssh_diffiehellman_group14;
		break;
	      case KEX_DHGROUP1:
		s->preferred_kex[s->n_preferred_kex++] =
		    &ssh_diffiehellman_group1;
		break;
	      case KEX_RSA:
		s->preferred_kex[s->n_preferred_kex++] =
		    &ssh_rsa_kex;
		break;
	      case KEX_WARN:
		/* Flag for later. Don't bother if it's the last in
		 * the list. */
		if (i < KEX_MAX - 1) {
		    s->preferred_kex[s->n_preferred_kex++] = NULL;
		}
		break;
	    }
	}

	/*
	 * Set up the preferred ciphers. (NULL => warn below here)
	 */
	s->n_preferred_ciphers = 0;
	for (i = 0; i < CIPHER_MAX; i++) {
	    switch (ssh->cfg.ssh_cipherlist[i]) {
	      case CIPHER_BLOWFISH:
		s->preferred_ciphers[s->n_preferred_ciphers++] = &ssh2_blowfish;
		break;
	      case CIPHER_DES:
		if (ssh->cfg.ssh2_des_cbc) {
		    s->preferred_ciphers[s->n_preferred_ciphers++] = &ssh2_des;
		}
		break;
	      case CIPHER_3DES:
		s->preferred_ciphers[s->n_preferred_ciphers++] = &ssh2_3des;
		break;
	      case CIPHER_AES:
		s->preferred_ciphers[s->n_preferred_ciphers++] = &ssh2_aes;
		break;
	      case CIPHER_ARCFOUR:
		s->preferred_ciphers[s->n_preferred_ciphers++] = &ssh2_arcfour;
		break;
	      case CIPHER_WARN:
		/* Flag for later. Don't bother if it's the last in
		 * the list. */
		if (i < CIPHER_MAX - 1) {
		    s->preferred_ciphers[s->n_preferred_ciphers++] = NULL;
		}
		break;
	    }
	}

	/*
	 * Set up preferred compression.
	 */
	if (ssh->cfg.compression)
	    s->preferred_comp = &ssh_zlib;
	else
	    s->preferred_comp = &ssh_comp_none;

	/*
	 * Enable queueing of outgoing auth- or connection-layer
	 * packets while we are in the middle of a key exchange.
	 */
	ssh->queueing = TRUE;

	/*
	 * Flag that KEX is in progress.
	 */
	ssh->kex_in_progress = TRUE;

	/*
	 * Construct and send our key exchange packet.
	 */
	s->pktout = ssh2_pkt_init(SSH2_MSG_KEXINIT);
	for (i = 0; i < 16; i++)
	    ssh2_pkt_addbyte(s->pktout, (unsigned char) random_byte());
	/* List key exchange algorithms. */
	ssh2_pkt_addstring_start(s->pktout);
	commalist_started = 0;
	for (i = 0; i < s->n_preferred_kex; i++) {
	    const struct ssh_kexes *k = s->preferred_kex[i];
	    if (!k) continue;	       /* warning flag */
	    for (j = 0; j < k->nkexes; j++) {
		if (commalist_started)
		    ssh2_pkt_addstring_str(s->pktout, ",");
		ssh2_pkt_addstring_str(s->pktout, k->list[j]->name);
		commalist_started = 1;
	    }
	}
	/* List server host key algorithms. */
	ssh2_pkt_addstring_start(s->pktout);
	for (i = 0; i < lenof(hostkey_algs); i++) {
	    ssh2_pkt_addstring_str(s->pktout, hostkey_algs[i]->name);
	    if (i < lenof(hostkey_algs) - 1)
		ssh2_pkt_addstring_str(s->pktout, ",");
	}
	/* List client->server encryption algorithms. */
	ssh2_pkt_addstring_start(s->pktout);
	commalist_started = 0;
	for (i = 0; i < s->n_preferred_ciphers; i++) {
	    const struct ssh2_ciphers *c = s->preferred_ciphers[i];
	    if (!c) continue;	       /* warning flag */
	    for (j = 0; j < c->nciphers; j++) {
		if (commalist_started)
		    ssh2_pkt_addstring_str(s->pktout, ",");
		ssh2_pkt_addstring_str(s->pktout, c->list[j]->name);
		commalist_started = 1;
	    }
	}
	/* List server->client encryption algorithms. */
	ssh2_pkt_addstring_start(s->pktout);
	commalist_started = 0;
	for (i = 0; i < s->n_preferred_ciphers; i++) {
	    const struct ssh2_ciphers *c = s->preferred_ciphers[i];
	    if (!c) continue; /* warning flag */
	    for (j = 0; j < c->nciphers; j++) {
		if (commalist_started)
		    ssh2_pkt_addstring_str(s->pktout, ",");
		ssh2_pkt_addstring_str(s->pktout, c->list[j]->name);
		commalist_started = 1;
	    }
	}
	/* List client->server MAC algorithms. */
	ssh2_pkt_addstring_start(s->pktout);
	for (i = 0; i < s->nmacs; i++) {
	    ssh2_pkt_addstring_str(s->pktout, s->maclist[i]->name);
	    if (i < s->nmacs - 1)
		ssh2_pkt_addstring_str(s->pktout, ",");
	}
	/* List server->client MAC algorithms. */
	ssh2_pkt_addstring_start(s->pktout);
	for (i = 0; i < s->nmacs; i++) {
	    ssh2_pkt_addstring_str(s->pktout, s->maclist[i]->name);
	    if (i < s->nmacs - 1)
		ssh2_pkt_addstring_str(s->pktout, ",");
	}
	/* List client->server compression algorithms,
	 * then server->client compression algorithms. (We use the
	 * same set twice.) */
	for (j = 0; j < 2; j++) {
	    ssh2_pkt_addstring_start(s->pktout);
	    assert(lenof(compressions) > 1);
	    /* Prefer non-delayed versions */
	    ssh2_pkt_addstring_str(s->pktout, s->preferred_comp->name);
	    /* We don't even list delayed versions of algorithms until
	     * they're allowed to be used, to avoid a race. See the end of
	     * this function. */
	    if (s->userauth_succeeded && s->preferred_comp->delayed_name) {
		ssh2_pkt_addstring_str(s->pktout, ",");
		ssh2_pkt_addstring_str(s->pktout,
				       s->preferred_comp->delayed_name);
	    }
	    for (i = 0; i < lenof(compressions); i++) {
		const struct ssh_compress *c = compressions[i];
		if (c != s->preferred_comp) {
		    ssh2_pkt_addstring_str(s->pktout, ",");
		    ssh2_pkt_addstring_str(s->pktout, c->name);
		    if (s->userauth_succeeded && c->delayed_name) {
			ssh2_pkt_addstring_str(s->pktout, ",");
			ssh2_pkt_addstring_str(s->pktout, c->delayed_name);
		    }
		}
	    }
	}
	/* List client->server languages. Empty list. */
	ssh2_pkt_addstring_start(s->pktout);
	/* List server->client languages. Empty list. */
	ssh2_pkt_addstring_start(s->pktout);
	/* First KEX packet does _not_ follow, because we're not that brave. */
	ssh2_pkt_addbool(s->pktout, FALSE);
	/* Reserved. */
	ssh2_pkt_adduint32(s->pktout, 0);
    }

    s->our_kexinitlen = s->pktout->length - 5;
    s->our_kexinit = snewn(s->our_kexinitlen, unsigned char);
    memcpy(s->our_kexinit, s->pktout->data + 5, s->our_kexinitlen); 

    ssh2_pkt_send_noqueue(ssh, s->pktout);

    if (!pktin)
	crWaitUntil(pktin);

    /*
     * Now examine the other side's KEXINIT to see what we're up
     * to.
     */
    {
	char *str, *preferred;
	int i, j, len;

	if (pktin->type != SSH2_MSG_KEXINIT) {
	    bombout(("expected key exchange packet from server"));
	    crStop(0);
	}
	ssh->kex = NULL;
	ssh->hostkey = NULL;
	s->cscipher_tobe = NULL;
	s->sccipher_tobe = NULL;
	s->csmac_tobe = NULL;
	s->scmac_tobe = NULL;
	s->cscomp_tobe = NULL;
	s->sccomp_tobe = NULL;
	s->warn_kex = s->warn_cscipher = s->warn_sccipher = FALSE;

	pktin->savedpos += 16;	        /* skip garbage cookie */
	ssh_pkt_getstring(pktin, &str, &len);    /* key exchange algorithms */

	preferred = NULL;
	for (i = 0; i < s->n_preferred_kex; i++) {
	    const struct ssh_kexes *k = s->preferred_kex[i];
	    if (!k) {
		s->warn_kex = TRUE;
	    } else {
		for (j = 0; j < k->nkexes; j++) {
		    if (!preferred) preferred = k->list[j]->name;
		    if (in_commasep_string(k->list[j]->name, str, len)) {
			ssh->kex = k->list[j];
			break;
		    }
		}
	    }
	    if (ssh->kex)
		break;
	}
	if (!ssh->kex) {
	    bombout(("Couldn't agree a key exchange algorithm (available: %s)",
		     str ? str : "(null)"));
	    crStop(0);
	}
	/*
	 * Note that the server's guess is considered wrong if it doesn't match
	 * the first algorithm in our list, even if it's still the algorithm
	 * we end up using.
	 */
	s->guessok = first_in_commasep_string(preferred, str, len);
	ssh_pkt_getstring(pktin, &str, &len);    /* host key algorithms */
	for (i = 0; i < lenof(hostkey_algs); i++) {
	    if (in_commasep_string(hostkey_algs[i]->name, str, len)) {
		ssh->hostkey = hostkey_algs[i];
		break;
	    }
	}
	s->guessok = s->guessok &&
	    first_in_commasep_string(hostkey_algs[0]->name, str, len);
	ssh_pkt_getstring(pktin, &str, &len);    /* client->server cipher */
	for (i = 0; i < s->n_preferred_ciphers; i++) {
	    const struct ssh2_ciphers *c = s->preferred_ciphers[i];
	    if (!c) {
		s->warn_cscipher = TRUE;
	    } else {
		for (j = 0; j < c->nciphers; j++) {
		    if (in_commasep_string(c->list[j]->name, str, len)) {
			s->cscipher_tobe = c->list[j];
			break;
		    }
		}
	    }
	    if (s->cscipher_tobe)
		break;
	}
	if (!s->cscipher_tobe) {
	    bombout(("Couldn't agree a client-to-server cipher (available: %s)",
		     str ? str : "(null)"));
	    crStop(0);
	}

	ssh_pkt_getstring(pktin, &str, &len);    /* server->client cipher */
	for (i = 0; i < s->n_preferred_ciphers; i++) {
	    const struct ssh2_ciphers *c = s->preferred_ciphers[i];
	    if (!c) {
		s->warn_sccipher = TRUE;
	    } else {
		for (j = 0; j < c->nciphers; j++) {
		    if (in_commasep_string(c->list[j]->name, str, len)) {
			s->sccipher_tobe = c->list[j];
			break;
		    }
		}
	    }
	    if (s->sccipher_tobe)
		break;
	}
	if (!s->sccipher_tobe) {
	    bombout(("Couldn't agree a server-to-client cipher (available: %s)",
		     str ? str : "(null)"));
	    crStop(0);
	}

	ssh_pkt_getstring(pktin, &str, &len);    /* client->server mac */
	for (i = 0; i < s->nmacs; i++) {
	    if (in_commasep_string(s->maclist[i]->name, str, len)) {
		s->csmac_tobe = s->maclist[i];
		break;
	    }
	}
	ssh_pkt_getstring(pktin, &str, &len);    /* server->client mac */
	for (i = 0; i < s->nmacs; i++) {
	    if (in_commasep_string(s->maclist[i]->name, str, len)) {
		s->scmac_tobe = s->maclist[i];
		break;
	    }
	}
	ssh_pkt_getstring(pktin, &str, &len);  /* client->server compression */
	for (i = 0; i < lenof(compressions) + 1; i++) {
	    const struct ssh_compress *c =
		i == 0 ? s->preferred_comp : compressions[i - 1];
	    if (in_commasep_string(c->name, str, len)) {
		s->cscomp_tobe = c;
		break;
	    } else if (in_commasep_string(c->delayed_name, str, len)) {
		if (s->userauth_succeeded) {
		    s->cscomp_tobe = c;
		    break;
		} else {
		    s->pending_compression = TRUE;  /* try this later */
		}
	    }
	}
	ssh_pkt_getstring(pktin, &str, &len);  /* server->client compression */
	for (i = 0; i < lenof(compressions) + 1; i++) {
	    const struct ssh_compress *c =
		i == 0 ? s->preferred_comp : compressions[i - 1];
	    if (in_commasep_string(c->name, str, len)) {
		s->sccomp_tobe = c;
		break;
	    } else if (in_commasep_string(c->delayed_name, str, len)) {
		if (s->userauth_succeeded) {
		    s->sccomp_tobe = c;
		    break;
		} else {
		    s->pending_compression = TRUE;  /* try this later */
		}
	    }
	}
	if (s->pending_compression) {
	    logevent("Server supports delayed compression; "
		     "will try this later");
	}
	ssh_pkt_getstring(pktin, &str, &len);  /* client->server language */
	ssh_pkt_getstring(pktin, &str, &len);  /* server->client language */
	s->ignorepkt = ssh2_pkt_getbool(pktin) && !s->guessok;

	if (s->warn_kex) {
	    ssh_set_frozen(ssh, 1);
	    s->dlgret = askalg(ssh->frontend, "key-exchange algorithm",
			       ssh->kex->name,
			       ssh_dialog_callback, ssh);
	    if (s->dlgret < 0) {
		do {
		    crReturn(0);
		    if (pktin) {
			bombout(("Unexpected data from server while"
				 " waiting for user response"));
			crStop(0);
		    }
		} while (pktin || inlen > 0);
		s->dlgret = ssh->user_response;
	    }
	    ssh_set_frozen(ssh, 0);
	    if (s->dlgret == 0) {
		ssh_disconnect(ssh, "User aborted at kex warning", NULL,
			       0, TRUE);
		crStop(0);
	    }
	}

	if (s->warn_cscipher) {
	    ssh_set_frozen(ssh, 1);
	    s->dlgret = askalg(ssh->frontend,
			       "client-to-server cipher",
			       s->cscipher_tobe->name,
			       ssh_dialog_callback, ssh);
	    if (s->dlgret < 0) {
		do {
		    crReturn(0);
		    if (pktin) {
			bombout(("Unexpected data from server while"
				 " waiting for user response"));
			crStop(0);
		    }
		} while (pktin || inlen > 0);
		s->dlgret = ssh->user_response;
	    }
	    ssh_set_frozen(ssh, 0);
	    if (s->dlgret == 0) {
		ssh_disconnect(ssh, "User aborted at cipher warning", NULL,
			       0, TRUE);
		crStop(0);
	    }
	}

	if (s->warn_sccipher) {
	    ssh_set_frozen(ssh, 1);
	    s->dlgret = askalg(ssh->frontend,
			       "server-to-client cipher",
			       s->sccipher_tobe->name,
			       ssh_dialog_callback, ssh);
	    if (s->dlgret < 0) {
		do {
		    crReturn(0);
		    if (pktin) {
			bombout(("Unexpected data from server while"
				 " waiting for user response"));
			crStop(0);
		    }
		} while (pktin || inlen > 0);
		s->dlgret = ssh->user_response;
	    }
	    ssh_set_frozen(ssh, 0);
	    if (s->dlgret == 0) {
		ssh_disconnect(ssh, "User aborted at cipher warning", NULL,
			       0, TRUE);
		crStop(0);
	    }
	}

	ssh->exhash = ssh->kex->hash->init();
	hash_string(ssh->kex->hash, ssh->exhash, ssh->v_c, strlen(ssh->v_c));
	hash_string(ssh->kex->hash, ssh->exhash, ssh->v_s, strlen(ssh->v_s));
	hash_string(ssh->kex->hash, ssh->exhash,
	    s->our_kexinit, s->our_kexinitlen);
	sfree(s->our_kexinit);
	if (pktin->length > 5)
	    hash_string(ssh->kex->hash, ssh->exhash,
		pktin->data + 5, pktin->length - 5);

	if (s->ignorepkt) /* first_kex_packet_follows */
	    crWaitUntil(pktin);                /* Ignore packet */
    }

    if (ssh->kex->main_type == KEXTYPE_DH) {
        /*
         * Work out the number of bits of key we will need from the
         * key exchange. We start with the maximum key length of
         * either cipher...
         */
        {
            int csbits, scbits;

            csbits = s->cscipher_tobe->keylen;
            scbits = s->sccipher_tobe->keylen;
            s->nbits = (csbits > scbits ? csbits : scbits);
        }
        /* The keys only have hlen-bit entropy, since they're based on
         * a hash. So cap the key size at hlen bits. */
        if (s->nbits > ssh->kex->hash->hlen * 8)
            s->nbits = ssh->kex->hash->hlen * 8;

        /*
         * If we're doing Diffie-Hellman group exchange, start by
         * requesting a group.
         */
        if (!ssh->kex->pdata) {
            logevent("Doing Diffie-Hellman group exchange");
            ssh->pkt_kctx = SSH2_PKTCTX_DHGEX;
            /*
             * Work out how big a DH group we will need to allow that
             * much data.
             */
            s->pbits = 512 << ((s->nbits - 1) / 64);
            s->pktout = ssh2_pkt_init(SSH2_MSG_KEX_DH_GEX_REQUEST);
            ssh2_pkt_adduint32(s->pktout, s->pbits);
            ssh2_pkt_send_noqueue(ssh, s->pktout);

            crWaitUntil(pktin);
            if (pktin->type != SSH2_MSG_KEX_DH_GEX_GROUP) {
                bombout(("expected key exchange group packet from server"));
                crStop(0);
            }
            s->p = ssh2_pkt_getmp(pktin);
            s->g = ssh2_pkt_getmp(pktin);
            if (!s->p || !s->g) {
                bombout(("unable to read mp-ints from incoming group packet"));
                crStop(0);
            }
            ssh->kex_ctx = dh_setup_gex(s->p, s->g);
            s->kex_init_value = SSH2_MSG_KEX_DH_GEX_INIT;
            s->kex_reply_value = SSH2_MSG_KEX_DH_GEX_REPLY;
        } else {
            ssh->pkt_kctx = SSH2_PKTCTX_DHGROUP;
            ssh->kex_ctx = dh_setup_group(ssh->kex);
            s->kex_init_value = SSH2_MSG_KEXDH_INIT;
            s->kex_reply_value = SSH2_MSG_KEXDH_REPLY;
            logeventf(ssh, "Using Diffie-Hellman with standard group \"%s\"",
                      ssh->kex->groupname);
        }

        logeventf(ssh, "Doing Diffie-Hellman key exchange with hash %s",
                  ssh->kex->hash->text_name);
        /*
         * Now generate and send e for Diffie-Hellman.
         */
        set_busy_status(ssh->frontend, BUSY_CPU); /* this can take a while */
        s->e = dh_create_e(ssh->kex_ctx, s->nbits * 2);
        s->pktout = ssh2_pkt_init(s->kex_init_value);
        ssh2_pkt_addmp(s->pktout, s->e);
        ssh2_pkt_send_noqueue(ssh, s->pktout);

        set_busy_status(ssh->frontend, BUSY_WAITING); /* wait for server */
        crWaitUntil(pktin);
        if (pktin->type != s->kex_reply_value) {
            bombout(("expected key exchange reply packet from server"));
            crStop(0);
        }
        set_busy_status(ssh->frontend, BUSY_CPU); /* cogitate */
        ssh_pkt_getstring(pktin, &s->hostkeydata, &s->hostkeylen);
        s->hkey = ssh->hostkey->newkey(s->hostkeydata, s->hostkeylen);
        s->f = ssh2_pkt_getmp(pktin);
        if (!s->f) {
            bombout(("unable to parse key exchange reply packet"));
            crStop(0);
        }
        ssh_pkt_getstring(pktin, &s->sigdata, &s->siglen);

        s->K = dh_find_K(ssh->kex_ctx, s->f);

        /* We assume everything from now on will be quick, and it might
         * involve user interaction. */
        set_busy_status(ssh->frontend, BUSY_NOT);

        hash_string(ssh->kex->hash, ssh->exhash, s->hostkeydata, s->hostkeylen);
        if (!ssh->kex->pdata) {
            hash_uint32(ssh->kex->hash, ssh->exhash, s->pbits);
            hash_mpint(ssh->kex->hash, ssh->exhash, s->p);
            hash_mpint(ssh->kex->hash, ssh->exhash, s->g);
        }
        hash_mpint(ssh->kex->hash, ssh->exhash, s->e);
        hash_mpint(ssh->kex->hash, ssh->exhash, s->f);

        dh_cleanup(ssh->kex_ctx);
        freebn(s->f);
        if (!ssh->kex->pdata) {
            freebn(s->g);
            freebn(s->p);
        }
    } else {
	logeventf(ssh, "Doing RSA key exchange with hash %s",
		  ssh->kex->hash->text_name);
	ssh->pkt_kctx = SSH2_PKTCTX_RSAKEX;
        /*
         * RSA key exchange. First expect a KEXRSA_PUBKEY packet
         * from the server.
         */
        crWaitUntil(pktin);
        if (pktin->type != SSH2_MSG_KEXRSA_PUBKEY) {
            bombout(("expected RSA public key packet from server"));
            crStop(0);
        }

        ssh_pkt_getstring(pktin, &s->hostkeydata, &s->hostkeylen);
        hash_string(ssh->kex->hash, ssh->exhash,
		    s->hostkeydata, s->hostkeylen);
	s->hkey = ssh->hostkey->newkey(s->hostkeydata, s->hostkeylen);

        {
            char *keydata;
            ssh_pkt_getstring(pktin, &keydata, &s->rsakeylen);
            s->rsakeydata = snewn(s->rsakeylen, char);
            memcpy(s->rsakeydata, keydata, s->rsakeylen);
        }

        s->rsakey = ssh_rsakex_newkey(s->rsakeydata, s->rsakeylen);
        if (!s->rsakey) {
            sfree(s->rsakeydata);
            bombout(("unable to parse RSA public key from server"));
            crStop(0);
        }

        hash_string(ssh->kex->hash, ssh->exhash, s->rsakeydata, s->rsakeylen);

        /*
         * Next, set up a shared secret K, of precisely KLEN -
         * 2*HLEN - 49 bits, where KLEN is the bit length of the
         * RSA key modulus and HLEN is the bit length of the hash
         * we're using.
         */
        {
            int klen = ssh_rsakex_klen(s->rsakey);
            int nbits = klen - (2*ssh->kex->hash->hlen*8 + 49);
            int i, byte = 0;
            unsigned char *kstr1, *kstr2, *outstr;
            int kstr1len, kstr2len, outstrlen;

            s->K = bn_power_2(nbits - 1);

            for (i = 0; i < nbits; i++) {
                if ((i & 7) == 0) {
                    byte = random_byte();
                }
                bignum_set_bit(s->K, i, (byte >> (i & 7)) & 1);
            }

            /*
             * Encode this as an mpint.
             */
            kstr1 = ssh2_mpint_fmt(s->K, &kstr1len);
            kstr2 = snewn(kstr2len = 4 + kstr1len, unsigned char);
            PUT_32BIT(kstr2, kstr1len);
            memcpy(kstr2 + 4, kstr1, kstr1len);

            /*
             * Encrypt it with the given RSA key.
             */
            outstrlen = (klen + 7) / 8;
            outstr = snewn(outstrlen, unsigned char);
            ssh_rsakex_encrypt(ssh->kex->hash, kstr2, kstr2len,
			       outstr, outstrlen, s->rsakey);

            /*
             * And send it off in a return packet.
             */
            s->pktout = ssh2_pkt_init(SSH2_MSG_KEXRSA_SECRET);
            ssh2_pkt_addstring_start(s->pktout);
            ssh2_pkt_addstring_data(s->pktout, (char *)outstr, outstrlen);
            ssh2_pkt_send_noqueue(ssh, s->pktout);

	    hash_string(ssh->kex->hash, ssh->exhash, outstr, outstrlen);

            sfree(kstr2);
            sfree(kstr1);
            sfree(outstr);
        }

        ssh_rsakex_freekey(s->rsakey);

        crWaitUntil(pktin);
        if (pktin->type != SSH2_MSG_KEXRSA_DONE) {
            sfree(s->rsakeydata);
            bombout(("expected signature packet from server"));
            crStop(0);
        }

        ssh_pkt_getstring(pktin, &s->sigdata, &s->siglen);

        sfree(s->rsakeydata);
    }

    hash_mpint(ssh->kex->hash, ssh->exhash, s->K);
    assert(ssh->kex->hash->hlen <= sizeof(s->exchange_hash));
    ssh->kex->hash->final(ssh->exhash, s->exchange_hash);

    ssh->kex_ctx = NULL;

#if 0
    debug(("Exchange hash is:\n"));
    dmemdump(s->exchange_hash, ssh->kex->hash->hlen);
#endif

    if (!s->hkey ||
	!ssh->hostkey->verifysig(s->hkey, s->sigdata, s->siglen,
				 (char *)s->exchange_hash,
				 ssh->kex->hash->hlen)) {
	bombout(("Server's host key did not match the signature supplied"));
	crStop(0);
    }

    /*
     * Authenticate remote host: verify host key. (We've already
     * checked the signature of the exchange hash.)
     */
    s->keystr = ssh->hostkey->fmtkey(s->hkey);
    s->fingerprint = ssh->hostkey->fingerprint(s->hkey);
    ssh_set_frozen(ssh, 1);
    s->dlgret = verify_ssh_host_key(ssh->frontend,
                                    ssh->savedhost, ssh->savedport,
                                    ssh->hostkey->keytype, s->keystr,
				    s->fingerprint,
                                    ssh_dialog_callback, ssh);
    if (s->dlgret < 0) {
        do {
            crReturn(0);
            if (pktin) {
                bombout(("Unexpected data from server while waiting"
                         " for user host key response"));
                    crStop(0);
            }
        } while (pktin || inlen > 0);
        s->dlgret = ssh->user_response;
    }
    ssh_set_frozen(ssh, 0);
    if (s->dlgret == 0) {
	ssh_disconnect(ssh, "User aborted at host key verification", NULL,
		       0, TRUE);
        crStop(0);
    }
    if (!s->got_session_id) {     /* don't bother logging this in rekeys */
	logevent("Host key fingerprint is:");
	logevent(s->fingerprint);
    }
    sfree(s->fingerprint);
    sfree(s->keystr);
    ssh->hostkey->freekey(s->hkey);

    /*
     * The exchange hash from the very first key exchange is also
     * the session id, used in session key construction and
     * authentication.
     */
    if (!s->got_session_id) {
	assert(sizeof(s->exchange_hash) <= sizeof(ssh->v2_session_id));
	memcpy(ssh->v2_session_id, s->exchange_hash,
	       sizeof(s->exchange_hash));
	ssh->v2_session_id_len = ssh->kex->hash->hlen;
	assert(ssh->v2_session_id_len <= sizeof(ssh->v2_session_id));
	s->got_session_id = TRUE;
    }

    /*
     * Send SSH2_MSG_NEWKEYS.
     */
    s->pktout = ssh2_pkt_init(SSH2_MSG_NEWKEYS);
    ssh2_pkt_send_noqueue(ssh, s->pktout);
    ssh->outgoing_data_size = 0;       /* start counting from here */

    /*
     * We've sent client NEWKEYS, so create and initialise
     * client-to-server session keys.
     */
    if (ssh->cs_cipher_ctx)
	ssh->cscipher->free_context(ssh->cs_cipher_ctx);
    ssh->cscipher = s->cscipher_tobe;
    ssh->cs_cipher_ctx = ssh->cscipher->make_context();

    if (ssh->cs_mac_ctx)
	ssh->csmac->free_context(ssh->cs_mac_ctx);
    ssh->csmac = s->csmac_tobe;
    ssh->cs_mac_ctx = ssh->csmac->make_context();

    if (ssh->cs_comp_ctx)
	ssh->cscomp->compress_cleanup(ssh->cs_comp_ctx);
    ssh->cscomp = s->cscomp_tobe;
    ssh->cs_comp_ctx = ssh->cscomp->compress_init();

    /*
     * Set IVs on client-to-server keys. Here we use the exchange
     * hash from the _first_ key exchange.
     */
    {
	unsigned char keyspace[SSH2_KEX_MAX_HASH_LEN * SSH2_MKKEY_ITERS];
	assert(sizeof(keyspace) >= ssh->kex->hash->hlen * SSH2_MKKEY_ITERS);
	ssh2_mkkey(ssh,s->K,s->exchange_hash,'C',keyspace);
	assert((ssh->cscipher->keylen+7) / 8 <=
	       ssh->kex->hash->hlen * SSH2_MKKEY_ITERS);
	ssh->cscipher->setkey(ssh->cs_cipher_ctx, keyspace);
	ssh2_mkkey(ssh,s->K,s->exchange_hash,'A',keyspace);
	assert(ssh->cscipher->blksize <=
	       ssh->kex->hash->hlen * SSH2_MKKEY_ITERS);
	ssh->cscipher->setiv(ssh->cs_cipher_ctx, keyspace);
	ssh2_mkkey(ssh,s->K,s->exchange_hash,'E',keyspace);
	assert(ssh->csmac->len <=
	       ssh->kex->hash->hlen * SSH2_MKKEY_ITERS);
	ssh->csmac->setkey(ssh->cs_mac_ctx, keyspace);
	memset(keyspace, 0, sizeof(keyspace));
    }

    logeventf(ssh, "Initialised %.200s client->server encryption",
	      ssh->cscipher->text_name);
    logeventf(ssh, "Initialised %.200s client->server MAC algorithm",
	      ssh->csmac->text_name);
    if (ssh->cscomp->text_name)
	logeventf(ssh, "Initialised %s compression",
		  ssh->cscomp->text_name);

    /*
     * Now our end of the key exchange is complete, we can send all
     * our queued higher-layer packets.
     */
    ssh->queueing = FALSE;
    ssh2_pkt_queuesend(ssh);

    /*
     * Expect SSH2_MSG_NEWKEYS from server.
     */
    crWaitUntil(pktin);
    if (pktin->type != SSH2_MSG_NEWKEYS) {
	bombout(("expected new-keys packet from server"));
	crStop(0);
    }
    ssh->incoming_data_size = 0;       /* start counting from here */

    /*
     * We've seen server NEWKEYS, so create and initialise
     * server-to-client session keys.
     */
    if (ssh->sc_cipher_ctx)
	ssh->sccipher->free_context(ssh->sc_cipher_ctx);
    ssh->sccipher = s->sccipher_tobe;
    ssh->sc_cipher_ctx = ssh->sccipher->make_context();

    if (ssh->sc_mac_ctx)
	ssh->scmac->free_context(ssh->sc_mac_ctx);
    ssh->scmac = s->scmac_tobe;
    ssh->sc_mac_ctx = ssh->scmac->make_context();

    if (ssh->sc_comp_ctx)
	ssh->sccomp->decompress_cleanup(ssh->sc_comp_ctx);
    ssh->sccomp = s->sccomp_tobe;
    ssh->sc_comp_ctx = ssh->sccomp->decompress_init();

    /*
     * Set IVs on server-to-client keys. Here we use the exchange
     * hash from the _first_ key exchange.
     */
    {
	unsigned char keyspace[SSH2_KEX_MAX_HASH_LEN * SSH2_MKKEY_ITERS];
	assert(sizeof(keyspace) >= ssh->kex->hash->hlen * SSH2_MKKEY_ITERS);
	ssh2_mkkey(ssh,s->K,s->exchange_hash,'D',keyspace);
	assert((ssh->sccipher->keylen+7) / 8 <=
	       ssh->kex->hash->hlen * SSH2_MKKEY_ITERS);
	ssh->sccipher->setkey(ssh->sc_cipher_ctx, keyspace);
	ssh2_mkkey(ssh,s->K,s->exchange_hash,'B',keyspace);
	assert(ssh->sccipher->blksize <=
	       ssh->kex->hash->hlen * SSH2_MKKEY_ITERS);
	ssh->sccipher->setiv(ssh->sc_cipher_ctx, keyspace);
	ssh2_mkkey(ssh,s->K,s->exchange_hash,'F',keyspace);
	assert(ssh->scmac->len <=
	       ssh->kex->hash->hlen * SSH2_MKKEY_ITERS);
	ssh->scmac->setkey(ssh->sc_mac_ctx, keyspace);
	memset(keyspace, 0, sizeof(keyspace));
    }
    logeventf(ssh, "Initialised %.200s server->client encryption",
	      ssh->sccipher->text_name);
    logeventf(ssh, "Initialised %.200s server->client MAC algorithm",
	      ssh->scmac->text_name);
    if (ssh->sccomp->text_name)
	logeventf(ssh, "Initialised %s decompression",
		  ssh->sccomp->text_name);

    /*
     * Free shared secret.
     */
    freebn(s->K);

    /*
     * Key exchange is over. Loop straight back round if we have a
     * deferred rekey reason.
     */
    if (ssh->deferred_rekey_reason) {
	logevent(ssh->deferred_rekey_reason);
	pktin = NULL;
	ssh->deferred_rekey_reason = NULL;
	goto begin_key_exchange;
    }

    /*
     * Otherwise, schedule a timer for our next rekey.
     */
    ssh->kex_in_progress = FALSE;
    ssh->last_rekey = GETTICKCOUNT();
    if (ssh->cfg.ssh_rekey_time != 0)
	ssh->next_rekey = schedule_timer(ssh->cfg.ssh_rekey_time*60*TICKSPERSEC,
					 ssh2_timer, ssh);

    /*
     * If this is the first key exchange phase, we must pass the
     * SSH2_MSG_NEWKEYS packet to the next layer, not because it
     * wants to see it but because it will need time to initialise
     * itself before it sees an actual packet. In subsequent key
     * exchange phases, we don't pass SSH2_MSG_NEWKEYS on, because
     * it would only confuse the layer above.
     */
    if (s->activated_authconn) {
	crReturn(0);
    }
    s->activated_authconn = TRUE;

    /*
     * Now we're encrypting. Begin returning 1 to the protocol main
     * function so that other things can run on top of the
     * transport. If we ever see a KEXINIT, we must go back to the
     * start.
     * 
     * We _also_ go back to the start if we see pktin==NULL and
     * inlen negative, because this is a special signal meaning
     * `initiate client-driven rekey', and `in' contains a message
     * giving the reason for the rekey.
     *
     * inlen==-1 means always initiate a rekey;
     * inlen==-2 means that userauth has completed successfully and
     *   we should consider rekeying (for delayed compression).
     */
    while (!((pktin && pktin->type == SSH2_MSG_KEXINIT) ||
	     (!pktin && inlen < 0))) {
        wait_for_rekey:
	crReturn(1);
    }
    if (pktin) {
	logevent("Server initiated key re-exchange");
    } else {
	if (inlen == -2) {
	    /* 
	     * authconn has seen a USERAUTH_SUCCEEDED. Time to enable
	     * delayed compression, if it's available.
	     *
	     * draft-miller-secsh-compression-delayed-00 says that you
	     * negotiate delayed compression in the first key exchange, and
	     * both sides start compressing when the server has sent
	     * USERAUTH_SUCCESS. This has a race condition -- the server
	     * can't know when the client has seen it, and thus which incoming
	     * packets it should treat as compressed.
	     *
	     * Instead, we do the initial key exchange without offering the
	     * delayed methods, but note if the server offers them; when we
	     * get here, if a delayed method was available that was higher
	     * on our list than what we got, we initiate a rekey in which we
	     * _do_ list the delayed methods (and hopefully get it as a
	     * result). Subsequent rekeys will do the same.
	     */
	    assert(!s->userauth_succeeded); /* should only happen once */
	    s->userauth_succeeded = TRUE;
	    if (!s->pending_compression)
		/* Can't see any point rekeying. */
		goto wait_for_rekey;       /* this is utterly horrid */
	    /* else fall through to rekey... */
	    s->pending_compression = FALSE;
	}
        /*
	 * Now we've decided to rekey.
	 *
         * Special case: if the server bug is set that doesn't
         * allow rekeying, we give a different log message and
         * continue waiting. (If such a server _initiates_ a rekey,
         * we process it anyway!)
         */
        if ((ssh->remote_bugs & BUG_SSH2_REKEY)) {
            logeventf(ssh, "Server bug prevents key re-exchange (%s)",
                      (char *)in);
            /* Reset the counters, so that at least this message doesn't
             * hit the event log _too_ often. */
            ssh->outgoing_data_size = 0;
            ssh->incoming_data_size = 0;
            if (ssh->cfg.ssh_rekey_time != 0) {
                ssh->next_rekey =
                    schedule_timer(ssh->cfg.ssh_rekey_time*60*TICKSPERSEC,
                                   ssh2_timer, ssh);
            }
            goto wait_for_rekey;       /* this is still utterly horrid */
        } else {
            logeventf(ssh, "Initiating key re-exchange (%s)", (char *)in);
        }
    }
    goto begin_key_exchange;

    crFinish(1);
}

/*
 * Add data to an SSH-2 channel output buffer.
 */
static void ssh2_add_channel_data(struct ssh_channel *c, char *buf,
				  int len)
{
    bufchain_add(&c->v.v2.outbuffer, buf, len);
}

/*
 * Attempt to send data on an SSH-2 channel.
 */
static int ssh2_try_send(struct ssh_channel *c)
{
    Ssh ssh = c->ssh;
    struct Packet *pktout;

    while (c->v.v2.remwindow > 0 && bufchain_size(&c->v.v2.outbuffer) > 0) {
	int len;
	void *data;
	bufchain_prefix(&c->v.v2.outbuffer, &data, &len);
	if ((unsigned)len > c->v.v2.remwindow)
	    len = c->v.v2.remwindow;
	if ((unsigned)len > c->v.v2.remmaxpkt)
	    len = c->v.v2.remmaxpkt;
	pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_DATA);
	ssh2_pkt_adduint32(pktout, c->remoteid);
	ssh2_pkt_addstring_start(pktout);
	dont_log_data(ssh, pktout, PKTLOG_OMIT);
	ssh2_pkt_addstring_data(pktout, data, len);
	end_log_omission(ssh, pktout);
	ssh2_pkt_send(ssh, pktout);
	bufchain_consume(&c->v.v2.outbuffer, len);
	c->v.v2.remwindow -= len;
    }

    /*
     * After having sent as much data as we can, return the amount
     * still buffered.
     */
    return bufchain_size(&c->v.v2.outbuffer);
}

static void ssh2_try_send_and_unthrottle(Ssh ssh, struct ssh_channel *c)
{
    int bufsize;
    if (c->closes)
	return;			       /* don't send on closing channels */
    bufsize = ssh2_try_send(c);
    if (bufsize == 0) {
	switch (c->type) {
	  case CHAN_MAINSESSION:
	    /* stdin need not receive an unthrottle
	     * notification since it will be polled */
	    break;
	  case CHAN_X11:
	    x11_unthrottle(c->u.x11.s);
	    break;
	  case CHAN_AGENT:
	    /* agent sockets are request/response and need no
	     * buffer management */
	    break;
	  case CHAN_SOCKDATA:
	    pfd_unthrottle(c->u.pfd.s);
	    break;
	}
    }

    /*
     * If we've emptied the channel's output buffer and there's a
     * pending close event, start the channel-closing procedure.
     */
    if (c->pending_close && bufchain_size(&c->v.v2.outbuffer) == 0) {
	struct Packet *pktout;
	pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_CLOSE);
	ssh2_pkt_adduint32(pktout, c->remoteid);
	ssh2_pkt_send(ssh, pktout);
	c->closes = 1;
	c->pending_close = FALSE;
    }
}

/*
 * Set up most of a new ssh_channel for SSH-2.
 */
static void ssh2_channel_init(struct ssh_channel *c)
{
    Ssh ssh = c->ssh;
    c->localid = alloc_channel_id(ssh);
    c->closes = 0;
    c->pending_close = FALSE;
    c->throttling_conn = FALSE;
    c->v.v2.locwindow = c->v.v2.locmaxwin = c->v.v2.remlocwin =
	ssh->cfg.ssh_simple ? OUR_V2_BIGWIN : OUR_V2_WINSIZE;
    c->v.v2.winadj_head = c->v.v2.winadj_tail = NULL;
    c->v.v2.throttle_state = UNTHROTTLED;
    bufchain_init(&c->v.v2.outbuffer);
}

/*
 * Potentially enlarge the window on an SSH-2 channel.
 */
static void ssh2_set_window(struct ssh_channel *c, int newwin)
{
    Ssh ssh = c->ssh;

    /*
     * Never send WINDOW_ADJUST for a channel that the remote side
     * already thinks it's closed; there's no point, since it won't
     * be sending any more data anyway.
     */
    if (c->closes != 0)
	return;

    /*
     * If the remote end has a habit of ignoring maxpkt, limit the
     * window so that it has no choice (assuming it doesn't ignore the
     * window as well).
     */
    if ((ssh->remote_bugs & BUG_SSH2_MAXPKT) && newwin > OUR_V2_MAXPKT)
	newwin = OUR_V2_MAXPKT;
	

    /*
     * Only send a WINDOW_ADJUST if there's significantly more window
     * available than the other end thinks there is.  This saves us
     * sending a WINDOW_ADJUST for every character in a shell session.
     *
     * "Significant" is arbitrarily defined as half the window size.
     */
    if (newwin / 2 >= c->v.v2.locwindow) {
	struct Packet *pktout;
	struct winadj *wa;

	/*
	 * In order to keep track of how much window the client
	 * actually has available, we'd like it to acknowledge each
	 * WINDOW_ADJUST.  We can't do that directly, so we accompany
	 * it with a CHANNEL_REQUEST that has to be acknowledged.
	 *
	 * This is only necessary if we're opening the window wide.
	 * If we're not, then throughput is being constrained by
	 * something other than the maximum window size anyway.
	 *
	 * We also only send this if the main channel has finished its
	 * initial CHANNEL_REQUESTs and installed the default
	 * CHANNEL_FAILURE handler, so as not to risk giving it
	 * unexpected CHANNEL_FAILUREs.
	 */
	if (newwin == c->v.v2.locmaxwin &&
	    ssh->packet_dispatch[SSH2_MSG_CHANNEL_FAILURE]) {
	    pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
	    ssh2_pkt_adduint32(pktout, c->remoteid);
	    ssh2_pkt_addstring(pktout, "winadj@putty.projects.tartarus.org");
	    ssh2_pkt_addbool(pktout, TRUE);
	    ssh2_pkt_send(ssh, pktout);

	    /*
	     * CHANNEL_FAILURE doesn't come with any indication of
	     * what message caused it, so we have to keep track of the
	     * outstanding CHANNEL_REQUESTs ourselves.
	     */
	    wa = snew(struct winadj);
	    wa->size = newwin - c->v.v2.locwindow;
	    wa->next = NULL;
	    if (!c->v.v2.winadj_head)
		c->v.v2.winadj_head = wa;
	    else
		c->v.v2.winadj_tail->next = wa;
	    c->v.v2.winadj_tail = wa;
	    if (c->v.v2.throttle_state != UNTHROTTLED)
		c->v.v2.throttle_state = UNTHROTTLING;
	} else {
	    /* Pretend the WINDOW_ADJUST was acked immediately. */
	    c->v.v2.remlocwin = newwin;
	    c->v.v2.throttle_state = THROTTLED;
	}
	pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_WINDOW_ADJUST);
	ssh2_pkt_adduint32(pktout, c->remoteid);
	ssh2_pkt_adduint32(pktout, newwin - c->v.v2.locwindow);
	ssh2_pkt_send(ssh, pktout);
	c->v.v2.locwindow = newwin;
    }
}

/*
 * Find the channel associated with a message.  If there's no channel,
 * or it's not properly open, make a noise about it and return NULL.
 */
static struct ssh_channel *ssh2_channel_msg(Ssh ssh, struct Packet *pktin)
{
    unsigned localid = ssh_pkt_getuint32(pktin);
    struct ssh_channel *c;

    c = find234(ssh->channels, &localid, ssh_channelfind);
    if (!c ||
	(c->halfopen && pktin->type != SSH2_MSG_CHANNEL_OPEN_CONFIRMATION &&
	 pktin->type != SSH2_MSG_CHANNEL_OPEN_FAILURE)) {
	char *buf = dupprintf("Received %s for %s channel %u",
			      ssh2_pkt_type(ssh->pkt_kctx, ssh->pkt_actx,
					    pktin->type),
			      c ? "half-open" : "nonexistent", localid);
	ssh_disconnect(ssh, NULL, buf, SSH2_DISCONNECT_PROTOCOL_ERROR, FALSE);
	sfree(buf);
	return NULL;
    }
    return c;
}

static int ssh2_handle_winadj_response(struct ssh_channel *c)
{
    struct winadj *wa = c->v.v2.winadj_head;
    if (!wa)
	return FALSE;
    c->v.v2.winadj_head = wa->next;
    c->v.v2.remlocwin += wa->size;
    sfree(wa);
    /*
     * winadj messages are only sent when the window is fully open, so
     * if we get an ack of one, we know any pending unthrottle is
     * complete.
     */
    if (c->v.v2.throttle_state == UNTHROTTLING)
	c->v.v2.throttle_state = UNTHROTTLED;
    return TRUE;
}

static void ssh2_msg_channel_success(Ssh ssh, struct Packet *pktin)
{
    /*
     * This should never get called.  All channel requests are either
     * sent with want_reply false, are sent before this handler gets
     * installed, or are "winadj@putty" requests, which servers should
     * never respond to with success.
     *
     * However, at least one server ("boks_sshd") is known to return
     * SUCCESS for channel requests it's never heard of, such as
     * "winadj@putty". Raised with foxt.com as bug 090916-090424, but
     * for the sake of a quiet life, we handle it just the same as the
     * expected FAILURE.
     */
    struct ssh_channel *c;

    c = ssh2_channel_msg(ssh, pktin);
    if (!c)
	return;
    if (!ssh2_handle_winadj_response(c))
	ssh_disconnect(ssh, NULL,
		       "Received unsolicited SSH_MSG_CHANNEL_SUCCESS",
		       SSH2_DISCONNECT_PROTOCOL_ERROR, FALSE);
}

static void ssh2_msg_channel_failure(Ssh ssh, struct Packet *pktin)
{
    /*
     * The only time this should get called is for "winadj@putty"
     * messages sent above.  All other channel requests are either
     * sent with want_reply false or are sent before this handler gets
     * installed.
     */
    struct ssh_channel *c;

    c = ssh2_channel_msg(ssh, pktin);
    if (!c)
	return;
    if (!ssh2_handle_winadj_response(c))
	ssh_disconnect(ssh, NULL,
		       "Received unsolicited SSH_MSG_CHANNEL_FAILURE",
		       SSH2_DISCONNECT_PROTOCOL_ERROR, FALSE);
}

static void ssh2_msg_channel_window_adjust(Ssh ssh, struct Packet *pktin)
{
    struct ssh_channel *c;
    c = ssh2_channel_msg(ssh, pktin);
    if (!c)
	return;
    if (!c->closes) {
	c->v.v2.remwindow += ssh_pkt_getuint32(pktin);
	ssh2_try_send_and_unthrottle(ssh, c);
    }
}

static void ssh2_msg_channel_data(Ssh ssh, struct Packet *pktin)
{
    char *data;
    int length;
    struct ssh_channel *c;
    c = ssh2_channel_msg(ssh, pktin);
    if (!c)
	return;
    if (pktin->type == SSH2_MSG_CHANNEL_EXTENDED_DATA &&
	ssh_pkt_getuint32(pktin) != SSH2_EXTENDED_DATA_STDERR)
	return;			       /* extended but not stderr */
    ssh_pkt_getstring(pktin, &data, &length);
    if (data) {
	int bufsize = 0;
	c->v.v2.locwindow -= length;
	c->v.v2.remlocwin -= length;
	switch (c->type) {
	  case CHAN_MAINSESSION:
	    bufsize =
		from_backend(ssh->frontend, pktin->type ==
			     SSH2_MSG_CHANNEL_EXTENDED_DATA,
			     data, length);
	    break;
	  case CHAN_X11:
	    bufsize = x11_send(c->u.x11.s, data, length);
	    break;
	  case CHAN_SOCKDATA:
	    bufsize = pfd_send(c->u.pfd.s, data, length);
	    break;
	  case CHAN_AGENT:
	    while (length > 0) {
		if (c->u.a.lensofar < 4) {
		    unsigned int l = min(4 - c->u.a.lensofar,
					 (unsigned)length);
		    memcpy(c->u.a.msglen + c->u.a.lensofar,
			   data, l);
		    data += l;
		    length -= l;
		    c->u.a.lensofar += l;
		}
		if (c->u.a.lensofar == 4) {
		    c->u.a.totallen =
			4 + GET_32BIT(c->u.a.msglen);
		    c->u.a.message = snewn(c->u.a.totallen,
					   unsigned char);
		    memcpy(c->u.a.message, c->u.a.msglen, 4);
		}
		if (c->u.a.lensofar >= 4 && length > 0) {
		    unsigned int l =
			min(c->u.a.totallen - c->u.a.lensofar,
			    (unsigned)length);
		    memcpy(c->u.a.message + c->u.a.lensofar,
			   data, l);
		    data += l;
		    length -= l;
		    c->u.a.lensofar += l;
		}
		if (c->u.a.lensofar == c->u.a.totallen) {
		    void *reply;
		    int replylen;
		    if (agent_query(c->u.a.message,
				    c->u.a.totallen,
				    &reply, &replylen,
				    ssh_agentf_callback, c))
			ssh_agentf_callback(c, reply, replylen);
		    sfree(c->u.a.message);
		    c->u.a.lensofar = 0;
		}
	    }
	    bufsize = 0;
	    break;
	}
	/*
	 * If it looks like the remote end hit the end of its window,
	 * and we didn't want it to do that, think about using a
	 * larger window.
	 */
	if (c->v.v2.remlocwin <= 0 && c->v.v2.throttle_state == UNTHROTTLED &&
	    c->v.v2.locmaxwin < 0x40000000)
	    c->v.v2.locmaxwin += OUR_V2_WINSIZE;
	/*
	 * If we are not buffering too much data,
	 * enlarge the window again at the remote side.
	 * If we are buffering too much, we may still
	 * need to adjust the window if the server's
	 * sent excess data.
	 */
	ssh2_set_window(c, bufsize < c->v.v2.locmaxwin ?
			c->v.v2.locmaxwin - bufsize : 0);
	/*
	 * If we're either buffering way too much data, or if we're
	 * buffering anything at all and we're in "simple" mode,
	 * throttle the whole channel.
	 */
	if ((bufsize > c->v.v2.locmaxwin ||
	     (ssh->cfg.ssh_simple && bufsize > 0)) &&
	    !c->throttling_conn) {
	    c->throttling_conn = 1;
	    ssh_throttle_conn(ssh, +1);
	}
    }
}

static void ssh2_msg_channel_eof(Ssh ssh, struct Packet *pktin)
{
    struct ssh_channel *c;

    c = ssh2_channel_msg(ssh, pktin);
    if (!c)
	return;

    if (c->type == CHAN_X11) {
	/*
	 * Remote EOF on an X11 channel means we should
	 * wrap up and close the channel ourselves.
	 */
	x11_close(c->u.x11.s);
	c->u.x11.s = NULL;
	sshfwd_close(c);
    } else if (c->type == CHAN_AGENT) {
	sshfwd_close(c);
    } else if (c->type == CHAN_SOCKDATA) {
	pfd_close(c->u.pfd.s);
	c->u.pfd.s = NULL;
	sshfwd_close(c);
    }
}

static void ssh2_msg_channel_close(Ssh ssh, struct Packet *pktin)
{
    struct ssh_channel *c;
    struct Packet *pktout;

    c = ssh2_channel_msg(ssh, pktin);
    if (!c)
	return;
    /* Do pre-close processing on the channel. */
    switch (c->type) {
      case CHAN_MAINSESSION:
	ssh->mainchan = NULL;
	update_specials_menu(ssh->frontend);
	break;
      case CHAN_X11:
	if (c->u.x11.s != NULL)
	    x11_close(c->u.x11.s);
	sshfwd_close(c);
	break;
      case CHAN_AGENT:
	sshfwd_close(c);
	break;
      case CHAN_SOCKDATA:
	if (c->u.pfd.s != NULL)
	    pfd_close(c->u.pfd.s);
	sshfwd_close(c);
	break;
    }
    if (c->closes == 0) {
	pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_CLOSE);
	ssh2_pkt_adduint32(pktout, c->remoteid);
	ssh2_pkt_send(ssh, pktout);
    }
    del234(ssh->channels, c);
    bufchain_clear(&c->v.v2.outbuffer);
    sfree(c);

    /*
     * See if that was the last channel left open.
     * (This is only our termination condition if we're
     * not running in -N mode.)
     */
    if (!ssh->cfg.ssh_no_shell && count234(ssh->channels) == 0) {
	/*
	 * We used to send SSH_MSG_DISCONNECT here,
	 * because I'd believed that _every_ conforming
	 * SSH-2 connection had to end with a disconnect
	 * being sent by at least one side; apparently
	 * I was wrong and it's perfectly OK to
	 * unceremoniously slam the connection shut
	 * when you're done, and indeed OpenSSH feels
	 * this is more polite than sending a
	 * DISCONNECT. So now we don't.
	 */
	ssh_disconnect(ssh, "All channels closed", NULL, 0, TRUE);
    }
}

static void ssh2_msg_channel_open_confirmation(Ssh ssh, struct Packet *pktin)
{
    struct ssh_channel *c;
    struct Packet *pktout;

    c = ssh2_channel_msg(ssh, pktin);
    if (!c)
	return;
    if (c->type != CHAN_SOCKDATA_DORMANT)
	return;			       /* dunno why they're confirming this */
    c->remoteid = ssh_pkt_getuint32(pktin);
    c->halfopen = FALSE;
    c->type = CHAN_SOCKDATA;
    c->v.v2.remwindow = ssh_pkt_getuint32(pktin);
    c->v.v2.remmaxpkt = ssh_pkt_getuint32(pktin);
    if (c->u.pfd.s)
	pfd_confirm(c->u.pfd.s);
    if (c->closes) {
	/*
	 * We have a pending close on this channel,
	 * which we decided on before the server acked
	 * the channel open. So now we know the
	 * remoteid, we can close it again.
	 */
	pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_CLOSE);
	ssh2_pkt_adduint32(pktout, c->remoteid);
	ssh2_pkt_send(ssh, pktout);
    }
}

static void ssh2_msg_channel_open_failure(Ssh ssh, struct Packet *pktin)
{
    static const char *const reasons[] = {
	"<unknown reason code>",
	    "Administratively prohibited",
	    "Connect failed",
	    "Unknown channel type",
	    "Resource shortage",
    };
    unsigned reason_code;
    char *reason_string;
    int reason_length;
    struct ssh_channel *c;
    c = ssh2_channel_msg(ssh, pktin);
    if (!c)
	return;
    if (c->type != CHAN_SOCKDATA_DORMANT)
	return;			       /* dunno why they're failing this */

    reason_code = ssh_pkt_getuint32(pktin);
    if (reason_code >= lenof(reasons))
	reason_code = 0; /* ensure reasons[reason_code] in range */
    ssh_pkt_getstring(pktin, &reason_string, &reason_length);
    logeventf(ssh, "Forwarded connection refused by server: %s [%.*s]",
	      reasons[reason_code], reason_length, reason_string);

    pfd_close(c->u.pfd.s);

    del234(ssh->channels, c);
    sfree(c);
}

static void ssh2_msg_channel_request(Ssh ssh, struct Packet *pktin)
{
    char *type;
    int typelen, want_reply;
    int reply = SSH2_MSG_CHANNEL_FAILURE; /* default */
    struct ssh_channel *c;
    struct Packet *pktout;

    c = ssh2_channel_msg(ssh, pktin);
    if (!c)
	return;
    ssh_pkt_getstring(pktin, &type, &typelen);
    want_reply = ssh2_pkt_getbool(pktin);

    /*
     * Having got the channel number, we now look at
     * the request type string to see if it's something
     * we recognise.
     */
    if (c == ssh->mainchan) {
	/*
	 * We recognise "exit-status" and "exit-signal" on
	 * the primary channel.
	 */
	if (typelen == 11 &&
	    !memcmp(type, "exit-status", 11)) {

	    ssh->exitcode = ssh_pkt_getuint32(pktin);
	    logeventf(ssh, "Server sent command exit status %d",
		      ssh->exitcode);
	    reply = SSH2_MSG_CHANNEL_SUCCESS;

	} else if (typelen == 11 &&
		   !memcmp(type, "exit-signal", 11)) {

	    int is_plausible = TRUE, is_int = FALSE;
	    char *fmt_sig = "", *fmt_msg = "";
	    char *msg;
	    int msglen = 0, core = FALSE;
	    /* ICK: older versions of OpenSSH (e.g. 3.4p1)
	     * provide an `int' for the signal, despite its
	     * having been a `string' in the drafts of RFC 4254 since at
	     * least 2001. (Fixed in session.c 1.147.) Try to
	     * infer which we can safely parse it as. */
	    {
		unsigned char *p = pktin->body +
		    pktin->savedpos;
		long len = pktin->length - pktin->savedpos;
		unsigned long num = GET_32BIT(p); /* what is it? */
		/* If it's 0, it hardly matters; assume string */
		if (num == 0) {
		    is_int = FALSE;
		} else {
		    int maybe_int = FALSE, maybe_str = FALSE;
#define CHECK_HYPOTHESIS(offset, result) \
    do { \
	long q = offset; \
	if (q >= 0 && q+4 <= len) { \
	    q = q + 4 + GET_32BIT(p+q); \
	    if (q >= 0 && q+4 <= len && \
		    ((q = q + 4 + GET_32BIT(p+q))!= 0) && q == len) \
		result = TRUE; \
	} \
    } while(0)
		    CHECK_HYPOTHESIS(4+1, maybe_int);
		    CHECK_HYPOTHESIS(4+num+1, maybe_str);
#undef CHECK_HYPOTHESIS
		    if (maybe_int && !maybe_str)
			is_int = TRUE;
		    else if (!maybe_int && maybe_str)
			is_int = FALSE;
		    else
			/* Crikey. Either or neither. Panic. */
			is_plausible = FALSE;
		}
	    }
	    ssh->exitcode = 128;       /* means `unknown signal' */
	    if (is_plausible) {
		if (is_int) {
		    /* Old non-standard OpenSSH. */
		    int signum = ssh_pkt_getuint32(pktin);
		    fmt_sig = dupprintf(" %d", signum);
		    ssh->exitcode = 128 + signum;
		} else {
		    /* As per RFC 4254. */
		    char *sig;
		    int siglen;
		    ssh_pkt_getstring(pktin, &sig, &siglen);
		    /* Signal name isn't supposed to be blank, but
		     * let's cope gracefully if it is. */
		    if (siglen) {
			fmt_sig = dupprintf(" \"%.*s\"",
					    siglen, sig);
		    }

		    /*
		     * Really hideous method of translating the
		     * signal description back into a locally
		     * meaningful number.
		     */

		    if (0)
			;
#define TRANSLATE_SIGNAL(s) \
    else if (siglen == lenof(#s)-1 && !memcmp(sig, #s, siglen)) \
        ssh->exitcode = 128 + SIG ## s
#ifdef SIGABRT
		    TRANSLATE_SIGNAL(ABRT);
#endif
#ifdef SIGALRM
		    TRANSLATE_SIGNAL(ALRM);
#endif
#ifdef SIGFPE
		    TRANSLATE_SIGNAL(FPE);
#endif
#ifdef SIGHUP
		    TRANSLATE_SIGNAL(HUP);
#endif
#ifdef SIGILL
		    TRANSLATE_SIGNAL(ILL);
#endif
#ifdef SIGINT
		    TRANSLATE_SIGNAL(INT);
#endif
#ifdef SIGKILL
		    TRANSLATE_SIGNAL(KILL);
#endif
#ifdef SIGPIPE
		    TRANSLATE_SIGNAL(PIPE);
#endif
#ifdef SIGQUIT
		    TRANSLATE_SIGNAL(QUIT);
#endif
#ifdef SIGSEGV
		    TRANSLATE_SIGNAL(SEGV);
#endif
#ifdef SIGTERM
		    TRANSLATE_SIGNAL(TERM);
#endif
#ifdef SIGUSR1
		    TRANSLATE_SIGNAL(USR1);
#endif
#ifdef SIGUSR2
		    TRANSLATE_SIGNAL(USR2);
#endif
#undef TRANSLATE_SIGNAL
		    else
			ssh->exitcode = 128;
		}
		core = ssh2_pkt_getbool(pktin);
		ssh_pkt_getstring(pktin, &msg, &msglen);
		if (msglen) {
		    fmt_msg = dupprintf(" (\"%.*s\")", msglen, msg);
		}
		/* ignore lang tag */
	    } /* else don't attempt to parse */
	    logeventf(ssh, "Server exited on signal%s%s%s",
		      fmt_sig, core ? " (core dumped)" : "",
		      fmt_msg);
	    if (*fmt_sig) sfree(fmt_sig);
	    if (*fmt_msg) sfree(fmt_msg);
	    reply = SSH2_MSG_CHANNEL_SUCCESS;

	}
    } else {
	/*
	 * This is a channel request we don't know
	 * about, so we now either ignore the request
	 * or respond with CHANNEL_FAILURE, depending
	 * on want_reply.
	 */
	reply = SSH2_MSG_CHANNEL_FAILURE;
    }
    if (want_reply) {
	pktout = ssh2_pkt_init(reply);
	ssh2_pkt_adduint32(pktout, c->remoteid);
	ssh2_pkt_send(ssh, pktout);
    }
}

static void ssh2_msg_global_request(Ssh ssh, struct Packet *pktin)
{
    char *type;
    int typelen, want_reply;
    struct Packet *pktout;

    ssh_pkt_getstring(pktin, &type, &typelen);
    want_reply = ssh2_pkt_getbool(pktin);

    /*
     * We currently don't support any global requests
     * at all, so we either ignore the request or
     * respond with REQUEST_FAILURE, depending on
     * want_reply.
     */
    if (want_reply) {
	pktout = ssh2_pkt_init(SSH2_MSG_REQUEST_FAILURE);
	ssh2_pkt_send(ssh, pktout);
    }
}

static void ssh2_msg_channel_open(Ssh ssh, struct Packet *pktin)
{
    char *type;
    int typelen;
    char *peeraddr;
    int peeraddrlen;
    int peerport;
    char *error = NULL;
    struct ssh_channel *c;
    unsigned remid, winsize, pktsize;
    struct Packet *pktout;

    ssh_pkt_getstring(pktin, &type, &typelen);
    c = snew(struct ssh_channel);
    c->ssh = ssh;

    remid = ssh_pkt_getuint32(pktin);
    winsize = ssh_pkt_getuint32(pktin);
    pktsize = ssh_pkt_getuint32(pktin);

    if (typelen == 3 && !memcmp(type, "x11", 3)) {
	char *addrstr;
	const char *x11err;

	ssh_pkt_getstring(pktin, &peeraddr, &peeraddrlen);
	addrstr = snewn(peeraddrlen+1, char);
	memcpy(addrstr, peeraddr, peeraddrlen);
	addrstr[peeraddrlen] = '\0';
	peerport = ssh_pkt_getuint32(pktin);

	logeventf(ssh, "Received X11 connect request from %s:%d",
		  addrstr, peerport);

	if (!ssh->X11_fwd_enabled)
	    error = "X11 forwarding is not enabled";
	else if ((x11err = x11_init(&c->u.x11.s, ssh->x11disp, c,
				    addrstr, peerport, &ssh->cfg)) != NULL) {
	    logeventf(ssh, "Local X11 connection failed: %s", x11err);
	    error = "Unable to open an X11 connection";
	} else {
	    logevent("Opening X11 forward connection succeeded");
	    c->type = CHAN_X11;
	}

	sfree(addrstr);
    } else if (typelen == 15 &&
	       !memcmp(type, "forwarded-tcpip", 15)) {
	struct ssh_rportfwd pf, *realpf;
	char *dummy;
	int dummylen;
	ssh_pkt_getstring(pktin, &dummy, &dummylen);/* skip address */
	pf.sport = ssh_pkt_getuint32(pktin);
	ssh_pkt_getstring(pktin, &peeraddr, &peeraddrlen);
	peerport = ssh_pkt_getuint32(pktin);
	realpf = find234(ssh->rportfwds, &pf, NULL);
	logeventf(ssh, "Received remote port %d open request "
		  "from %s:%d", pf.sport, peeraddr, peerport);
	if (realpf == NULL) {
	    error = "Remote port is not recognised";
	} else {
	    const char *e = pfd_newconnect(&c->u.pfd.s,
					   realpf->dhost,
					   realpf->dport, c,
					   &ssh->cfg,
					   realpf->pfrec->addressfamily);
	    logeventf(ssh, "Attempting to forward remote port to "
		      "%s:%d", realpf->dhost, realpf->dport);
	    if (e != NULL) {
		logeventf(ssh, "Port open failed: %s", e);
		error = "Port open failed";
	    } else {
		logevent("Forwarded port opened successfully");
		c->type = CHAN_SOCKDATA;
	    }
	}
    } else if (typelen == 22 &&
	       !memcmp(type, "auth-agent@openssh.com", 22)) {
	if (!ssh->agentfwd_enabled)
	    error = "Agent forwarding is not enabled";
	else {
	    c->type = CHAN_AGENT;	/* identify channel type */
	    c->u.a.lensofar = 0;
	}
    } else {
	error = "Unsupported channel type requested";
    }

    c->remoteid = remid;
    c->halfopen = FALSE;
    if (error) {
	pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_OPEN_FAILURE);
	ssh2_pkt_adduint32(pktout, c->remoteid);
	ssh2_pkt_adduint32(pktout, SSH2_OPEN_CONNECT_FAILED);
	ssh2_pkt_addstring(pktout, error);
	ssh2_pkt_addstring(pktout, "en");	/* language tag */
	ssh2_pkt_send(ssh, pktout);
	logeventf(ssh, "Rejected channel open: %s", error);
	sfree(c);
    } else {
	ssh2_channel_init(c);
	c->v.v2.remwindow = winsize;
	c->v.v2.remmaxpkt = pktsize;
	add234(ssh->channels, c);
	pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_OPEN_CONFIRMATION);
	ssh2_pkt_adduint32(pktout, c->remoteid);
	ssh2_pkt_adduint32(pktout, c->localid);
	ssh2_pkt_adduint32(pktout, c->v.v2.locwindow);
	ssh2_pkt_adduint32(pktout, OUR_V2_MAXPKT);	/* our max pkt size */
	ssh2_pkt_send(ssh, pktout);
    }
}

/*
 * Buffer banner messages for later display at some convenient point,
 * if we're going to display them.
 */
static void ssh2_msg_userauth_banner(Ssh ssh, struct Packet *pktin)
{
    /* Arbitrary limit to prevent unbounded inflation of buffer */
    if (ssh->cfg.ssh_show_banner &&
	bufchain_size(&ssh->banner) <= 131072) {
	char *banner = NULL;
	int size = 0;
	ssh_pkt_getstring(pktin, &banner, &size);
	if (banner)
	    bufchain_add(&ssh->banner, banner, size);
    }
}

/* Helper function to deal with sending tty modes for "pty-req" */
static void ssh2_send_ttymode(void *data, char *mode, char *val)
{
    struct Packet *pktout = (struct Packet *)data;
    int i = 0;
    unsigned int arg = 0;
    while (strcmp(mode, ssh_ttymodes[i].mode) != 0) i++;
    if (i == lenof(ssh_ttymodes)) return;
    switch (ssh_ttymodes[i].type) {
      case TTY_OP_CHAR:
	arg = ssh_tty_parse_specchar(val);
	break;
      case TTY_OP_BOOL:
	arg = ssh_tty_parse_boolean(val);
	break;
    }
    ssh2_pkt_addbyte(pktout, ssh_ttymodes[i].opcode);
    ssh2_pkt_adduint32(pktout, arg);
}

/*
 * Handle the SSH-2 userauth and connection layers.
 */
static void do_ssh2_authconn(Ssh ssh, unsigned char *in, int inlen,
			     struct Packet *pktin)
{
    struct do_ssh2_authconn_state {
	enum {
	    AUTH_TYPE_NONE,
		AUTH_TYPE_PUBLICKEY,
		AUTH_TYPE_PUBLICKEY_OFFER_LOUD,
		AUTH_TYPE_PUBLICKEY_OFFER_QUIET,
		AUTH_TYPE_PASSWORD,
	        AUTH_TYPE_GSSAPI,      /* always QUIET */
		AUTH_TYPE_KEYBOARD_INTERACTIVE,
		AUTH_TYPE_KEYBOARD_INTERACTIVE_QUIET
	} type;
	int done_service_req;
	int gotit, need_pw, can_pubkey, can_passwd, can_keyb_inter;
	int tried_pubkey_config, done_agent;
#ifndef NO_GSSAPI
	int can_gssapi;
	int tried_gssapi;
#endif
	int kbd_inter_refused;
	int we_are_in, userauth_success;
	prompts_t *cur_prompt;
	int num_prompts;
	char username[100];
	char *password;
	int got_username;
	void *publickey_blob;
	int publickey_bloblen;
	int publickey_encrypted;
	char *publickey_algorithm;
	char *publickey_comment;
	unsigned char agent_request[5], *agent_response, *agentp;
	int agent_responselen;
	unsigned char *pkblob_in_agent;
	int keyi, nkeys;
	char *pkblob, *alg, *commentp;
	int pklen, alglen, commentlen;
	int siglen, retlen, len;
	char *q, *agentreq, *ret;
	int try_send;
	int num_env, env_left, env_ok;
	struct Packet *pktout;
#ifndef NO_GSSAPI
	struct ssh_gss_library *gsslib;
	Ssh_gss_ctx gss_ctx;
	Ssh_gss_buf gss_buf;
	Ssh_gss_buf gss_rcvtok, gss_sndtok;
	Ssh_gss_name gss_srv_name;
	Ssh_gss_stat gss_stat;
#endif
    };
    crState(do_ssh2_authconn_state);

    crBegin(ssh->do_ssh2_authconn_crstate);

    s->done_service_req = FALSE;
    s->we_are_in = s->userauth_success = FALSE;
#ifndef NO_GSSAPI
    s->tried_gssapi = FALSE;
#endif

    if (!ssh->cfg.ssh_no_userauth) {
	/*
	 * Request userauth protocol, and await a response to it.
	 */
	s->pktout = ssh2_pkt_init(SSH2_MSG_SERVICE_REQUEST);
	ssh2_pkt_addstring(s->pktout, "ssh-userauth");
	ssh2_pkt_send(ssh, s->pktout);
	crWaitUntilV(pktin);
	if (pktin->type == SSH2_MSG_SERVICE_ACCEPT)
	    s->done_service_req = TRUE;
    }
    if (!s->done_service_req) {
	/*
	 * Request connection protocol directly, without authentication.
	 */
	s->pktout = ssh2_pkt_init(SSH2_MSG_SERVICE_REQUEST);
	ssh2_pkt_addstring(s->pktout, "ssh-connection");
	ssh2_pkt_send(ssh, s->pktout);
	crWaitUntilV(pktin);
	if (pktin->type == SSH2_MSG_SERVICE_ACCEPT) {
	    s->we_are_in = TRUE; /* no auth required */
	} else {
	    bombout(("Server refused service request"));
	    crStopV;
	}
    }

    /* Arrange to be able to deal with any BANNERs that come in.
     * (We do this now as packets may come in during the next bit.) */
    bufchain_init(&ssh->banner);
    ssh->packet_dispatch[SSH2_MSG_USERAUTH_BANNER] =
	ssh2_msg_userauth_banner;

    /*
     * Misc one-time setup for authentication.
     */
    s->publickey_blob = NULL;
    if (!s->we_are_in) {

	/*
	 * Load the public half of any configured public key file
	 * for later use.
	 */
	if (!filename_is_null(ssh->cfg.keyfile)) {
	    int keytype;
	    logeventf(ssh, "Reading private key file \"%.150s\"",
		      filename_to_str(&ssh->cfg.keyfile));
	    keytype = key_type(&ssh->cfg.keyfile);
	    if (keytype == SSH_KEYTYPE_SSH2) {
		const char *error;
		s->publickey_blob =
		    ssh2_userkey_loadpub(&ssh->cfg.keyfile,
					 &s->publickey_algorithm,
					 &s->publickey_bloblen, 
					 &s->publickey_comment, &error);
		if (s->publickey_blob) {
		    s->publickey_encrypted =
			ssh2_userkey_encrypted(&ssh->cfg.keyfile, NULL);
		} else {
		    char *msgbuf;
		    logeventf(ssh, "Unable to load private key (%s)", 
			      error);
		    msgbuf = dupprintf("Unable to load private key file "
				       "\"%.150s\" (%s)\r\n",
				       filename_to_str(&ssh->cfg.keyfile),
				       error);
		    c_write_str(ssh, msgbuf);
		    sfree(msgbuf);
		}
	    } else {
		char *msgbuf;
		logeventf(ssh, "Unable to use this key file (%s)",
			  key_type_to_str(keytype));
		msgbuf = dupprintf("Unable to use key file \"%.150s\""
				   " (%s)\r\n",
				   filename_to_str(&ssh->cfg.keyfile),
				   key_type_to_str(keytype));
		c_write_str(ssh, msgbuf);
		sfree(msgbuf);
		s->publickey_blob = NULL;
	    }
	}

	/*
	 * Find out about any keys Pageant has (but if there's a
	 * public key configured, filter out all others).
	 */
	s->nkeys = 0;
	s->agent_response = NULL;
	s->pkblob_in_agent = NULL;
	if (ssh->cfg.tryagent && agent_exists()) {

	    void *r;

	    logevent("Pageant is running. Requesting keys.");

	    /* Request the keys held by the agent. */
	    PUT_32BIT(s->agent_request, 1);
	    s->agent_request[4] = SSH2_AGENTC_REQUEST_IDENTITIES;
	    if (!agent_query(s->agent_request, 5, &r, &s->agent_responselen,
			     ssh_agent_callback, ssh)) {
		do {
		    crReturnV;
		    if (pktin) {
			bombout(("Unexpected data from server while"
				 " waiting for agent response"));
			crStopV;
		    }
		} while (pktin || inlen > 0);
		r = ssh->agent_response;
		s->agent_responselen = ssh->agent_response_len;
	    }
	    s->agent_response = (unsigned char *) r;
	    if (s->agent_response && s->agent_responselen >= 5 &&
		s->agent_response[4] == SSH2_AGENT_IDENTITIES_ANSWER) {
		int keyi;
		unsigned char *p;
		p = s->agent_response + 5;
		s->nkeys = GET_32BIT(p);
		p += 4;
		logeventf(ssh, "Pageant has %d SSH-2 keys", s->nkeys);
		if (s->publickey_blob) {
		    /* See if configured key is in agent. */
		    for (keyi = 0; keyi < s->nkeys; keyi++) {
			s->pklen = GET_32BIT(p);
			if (s->pklen == s->publickey_bloblen &&
			    !memcmp(p+4, s->publickey_blob,
				    s->publickey_bloblen)) {
			    logeventf(ssh, "Pageant key #%d matches "
				      "configured key file", keyi);
			    s->keyi = keyi;
			    s->pkblob_in_agent = p;
			    break;
			}
			p += 4 + s->pklen;
			p += GET_32BIT(p) + 4; /* comment */
		    }
		    if (!s->pkblob_in_agent) {
			logevent("Configured key file not in Pageant");
			s->nkeys = 0;
		    }
		}
	    } else {
                logevent("Failed to get reply from Pageant");
	    }
	}

    }

    /*
     * We repeat this whole loop, including the username prompt,
     * until we manage a successful authentication. If the user
     * types the wrong _password_, they can be sent back to the
     * beginning to try another username, if this is configured on.
     * (If they specify a username in the config, they are never
     * asked, even if they do give a wrong password.)
     * 
     * I think this best serves the needs of
     * 
     *  - the people who have no configuration, no keys, and just
     *    want to try repeated (username,password) pairs until they
     *    type both correctly
     * 
     *  - people who have keys and configuration but occasionally
     *    need to fall back to passwords
     * 
     *  - people with a key held in Pageant, who might not have
     *    logged in to a particular machine before; so they want to
     *    type a username, and then _either_ their key will be
     *    accepted, _or_ they will type a password. If they mistype
     *    the username they will want to be able to get back and
     *    retype it!
     */
    s->username[0] = '\0';
    s->got_username = FALSE;
    while (!s->we_are_in) {
	/*
	 * Get a username.
	 */
	if (s->got_username && !ssh->cfg.change_username) {
	    /*
	     * We got a username last time round this loop, and
	     * with change_username turned off we don't try to get
	     * it again.
	     */
	} else if (!get_remote_username(&ssh->cfg, s->username,
					sizeof(s->username))) {
	    int ret; /* need not be kept over crReturn */
	    s->cur_prompt = new_prompts(ssh->frontend);
	    s->cur_prompt->to_server = TRUE;
	    s->cur_prompt->name = dupstr("SSH login name");
	    add_prompt(s->cur_prompt, dupstr("login as: "), TRUE,
		       lenof(s->username)); 
	    ret = get_userpass_input(s->cur_prompt, NULL, 0);
	    while (ret < 0) {
		ssh->send_ok = 1;
		crWaitUntilV(!pktin);
		ret = get_userpass_input(s->cur_prompt, in, inlen);
		ssh->send_ok = 0;
	    }
	    if (!ret) {
		/*
		 * get_userpass_input() failed to get a username.
		 * Terminate.
		 */
		free_prompts(s->cur_prompt);
		ssh_disconnect(ssh, "No username provided", NULL, 0, TRUE);
		crStopV;
	    }
	    memcpy(s->username, s->cur_prompt->prompts[0]->result,
		   lenof(s->username));
	    free_prompts(s->cur_prompt);
	} else {
	    char *stuff;
	    if ((flags & FLAG_VERBOSE) || (flags & FLAG_INTERACTIVE)) {
		stuff = dupprintf("Using username \"%s\".\r\n", s->username);
		c_write_str(ssh, stuff);
		sfree(stuff);
	    }
	}
	s->got_username = TRUE;

	/*
	 * Send an authentication request using method "none": (a)
	 * just in case it succeeds, and (b) so that we know what
	 * authentication methods we can usefully try next.
	 */
	ssh->pkt_actx = SSH2_PKTCTX_NOAUTH;

	s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_REQUEST);
	ssh2_pkt_addstring(s->pktout, s->username);
	ssh2_pkt_addstring(s->pktout, "ssh-connection");/* service requested */
	ssh2_pkt_addstring(s->pktout, "none");    /* method */
	ssh2_pkt_send(ssh, s->pktout);
	s->type = AUTH_TYPE_NONE;
	s->gotit = FALSE;
	s->we_are_in = FALSE;

	s->tried_pubkey_config = FALSE;
	s->kbd_inter_refused = FALSE;

	/* Reset agent request state. */
	s->done_agent = FALSE;
	if (s->agent_response) {
	    if (s->pkblob_in_agent) {
		s->agentp = s->pkblob_in_agent;
	    } else {
		s->agentp = s->agent_response + 5 + 4;
		s->keyi = 0;
	    }
	}

	while (1) {
	    char *methods = NULL;
	    int methlen = 0;

	    /*
	     * Wait for the result of the last authentication request.
	     */
	    if (!s->gotit)
		crWaitUntilV(pktin);
	    /*
	     * Now is a convenient point to spew any banner material
	     * that we've accumulated. (This should ensure that when
	     * we exit the auth loop, we haven't any left to deal
	     * with.)
	     */
	    {
		int size = bufchain_size(&ssh->banner);
		/*
		 * Don't show the banner if we're operating in
		 * non-verbose non-interactive mode. (It's probably
		 * a script, which means nobody will read the
		 * banner _anyway_, and moreover the printing of
		 * the banner will screw up processing on the
		 * output of (say) plink.)
		 */
		if (size && (flags & (FLAG_VERBOSE | FLAG_INTERACTIVE))) {
		    char *banner = snewn(size, char);
		    bufchain_fetch(&ssh->banner, banner, size);
		    c_write_untrusted(ssh, banner, size);
		    sfree(banner);
		}
		bufchain_clear(&ssh->banner);
	    }
	    if (pktin->type == SSH2_MSG_USERAUTH_SUCCESS) {
		logevent("Access granted");
		s->we_are_in = s->userauth_success = TRUE;
		break;
	    }

	    if (pktin->type != SSH2_MSG_USERAUTH_FAILURE && s->type != AUTH_TYPE_GSSAPI) {
		bombout(("Strange packet received during authentication: "
			 "type %d", pktin->type));
		crStopV;
	    }

	    s->gotit = FALSE;

	    /*
	     * OK, we're now sitting on a USERAUTH_FAILURE message, so
	     * we can look at the string in it and know what we can
	     * helpfully try next.
	     */
	    if (pktin->type == SSH2_MSG_USERAUTH_FAILURE) {
		ssh_pkt_getstring(pktin, &methods, &methlen);
		if (!ssh2_pkt_getbool(pktin)) {
		    /*
		     * We have received an unequivocal Access
		     * Denied. This can translate to a variety of
		     * messages, or no message at all.
                     *
                     * For forms of authentication which are attempted
                     * implicitly, by which I mean without printing
                     * anything in the window indicating that we're
                     * trying them, we should never print 'Access
                     * denied'.
                     *
                     * If we do print a message saying that we're
                     * attempting some kind of authentication, it's OK
                     * to print a followup message saying it failed -
                     * but the message may sometimes be more specific
                     * than simply 'Access denied'.
                     *
		     * Additionally, if we'd just tried password
		     * authentication, we should break out of this
		     * whole loop so as to go back to the username
		     * prompt (iff we're configured to allow
		     * username change attempts).
		     */
		    if (s->type == AUTH_TYPE_NONE) {
			/* do nothing */
		    } else if (s->type == AUTH_TYPE_PUBLICKEY_OFFER_LOUD ||
			       s->type == AUTH_TYPE_PUBLICKEY_OFFER_QUIET) {
			if (s->type == AUTH_TYPE_PUBLICKEY_OFFER_LOUD)
			    c_write_str(ssh, "Server refused our key\r\n");
			logevent("Server refused our key");
                    } else if (s->type == AUTH_TYPE_PUBLICKEY) {
                        /* This _shouldn't_ happen except by a
                         * protocol bug causing client and server to
                         * disagree on what is a correct signature. */
                        c_write_str(ssh, "Server refused public-key signature"
                                    " despite accepting key!\r\n");
                        logevent("Server refused public-key signature"
                                 " despite accepting key!");
		    } else if (s->type==AUTH_TYPE_KEYBOARD_INTERACTIVE_QUIET) {
                        /* quiet, so no c_write */
                        logevent("Server refused keyboard-interactive authentication");
		    } else if (s->type==AUTH_TYPE_GSSAPI) {
			/* always quiet, so no c_write */
                        /* also, the code down in the GSSAPI block has
                         * already logged this in the Event Log */
		    } else if (s->type == AUTH_TYPE_KEYBOARD_INTERACTIVE) {
                        logevent("Keyboard-interactive authentication failed");
			c_write_str(ssh, "Access denied\r\n");
                    } else {
                        assert(s->type == AUTH_TYPE_PASSWORD);
                        logevent("Password authentication failed");
			c_write_str(ssh, "Access denied\r\n");

			if (ssh->cfg.change_username) {
			    /* XXX perhaps we should allow
			     * keyboard-interactive to do this too? */
			    s->we_are_in = FALSE;
			    break;
			}
		    }
		} else {
		    c_write_str(ssh, "Further authentication required\r\n");
		    logevent("Further authentication required");
		}

		s->can_pubkey =
		    in_commasep_string("publickey", methods, methlen);
		s->can_passwd =
		    in_commasep_string("password", methods, methlen);
		s->can_keyb_inter = ssh->cfg.try_ki_auth &&
		    in_commasep_string("keyboard-interactive", methods, methlen);
#ifndef NO_GSSAPI
		if (!ssh->gsslibs)
		    ssh->gsslibs = ssh_gss_setup(&ssh->cfg);
		s->can_gssapi = ssh->cfg.try_gssapi_auth &&
		    in_commasep_string("gssapi-with-mic", methods, methlen) &&
		    ssh->gsslibs->nlibraries > 0;
#endif
	    }

	    ssh->pkt_actx = SSH2_PKTCTX_NOAUTH;

	    if (s->can_pubkey && !s->done_agent && s->nkeys) {

		/*
		 * Attempt public-key authentication using a key from Pageant.
		 */

		ssh->pkt_actx = SSH2_PKTCTX_PUBLICKEY;

		logeventf(ssh, "Trying Pageant key #%d", s->keyi);

		/* Unpack key from agent response */
		s->pklen = GET_32BIT(s->agentp);
		s->agentp += 4;
		s->pkblob = (char *)s->agentp;
		s->agentp += s->pklen;
		s->alglen = GET_32BIT(s->pkblob);
		s->alg = s->pkblob + 4;
		s->commentlen = GET_32BIT(s->agentp);
		s->agentp += 4;
		s->commentp = (char *)s->agentp;
		s->agentp += s->commentlen;
		/* s->agentp now points at next key, if any */

		/* See if server will accept it */
		s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_REQUEST);
		ssh2_pkt_addstring(s->pktout, s->username);
		ssh2_pkt_addstring(s->pktout, "ssh-connection");
						    /* service requested */
		ssh2_pkt_addstring(s->pktout, "publickey");
						    /* method */
		ssh2_pkt_addbool(s->pktout, FALSE); /* no signature included */
		ssh2_pkt_addstring_start(s->pktout);
		ssh2_pkt_addstring_data(s->pktout, s->alg, s->alglen);
		ssh2_pkt_addstring_start(s->pktout);
		ssh2_pkt_addstring_data(s->pktout, s->pkblob, s->pklen);
		ssh2_pkt_send(ssh, s->pktout);
		s->type = AUTH_TYPE_PUBLICKEY_OFFER_QUIET;

		crWaitUntilV(pktin);
		if (pktin->type != SSH2_MSG_USERAUTH_PK_OK) {

		    /* Offer of key refused. */
		    s->gotit = TRUE;

		} else {
		    
		    void *vret;

		    if (flags & FLAG_VERBOSE) {
			c_write_str(ssh, "Authenticating with "
				    "public key \"");
			c_write(ssh, s->commentp, s->commentlen);
			c_write_str(ssh, "\" from agent\r\n");
		    }

		    /*
		     * Server is willing to accept the key.
		     * Construct a SIGN_REQUEST.
		     */
		    s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_REQUEST);
		    ssh2_pkt_addstring(s->pktout, s->username);
		    ssh2_pkt_addstring(s->pktout, "ssh-connection");
							/* service requested */
		    ssh2_pkt_addstring(s->pktout, "publickey");
							/* method */
		    ssh2_pkt_addbool(s->pktout, TRUE);  /* signature included */
		    ssh2_pkt_addstring_start(s->pktout);
		    ssh2_pkt_addstring_data(s->pktout, s->alg, s->alglen);
		    ssh2_pkt_addstring_start(s->pktout);
		    ssh2_pkt_addstring_data(s->pktout, s->pkblob, s->pklen);

		    /* Ask agent for signature. */
		    s->siglen = s->pktout->length - 5 + 4 +
			ssh->v2_session_id_len;
		    if (ssh->remote_bugs & BUG_SSH2_PK_SESSIONID)
			s->siglen -= 4;
		    s->len = 1;       /* message type */
		    s->len += 4 + s->pklen;	/* key blob */
		    s->len += 4 + s->siglen;	/* data to sign */
		    s->len += 4;      /* flags */
		    s->agentreq = snewn(4 + s->len, char);
		    PUT_32BIT(s->agentreq, s->len);
		    s->q = s->agentreq + 4;
		    *s->q++ = SSH2_AGENTC_SIGN_REQUEST;
		    PUT_32BIT(s->q, s->pklen);
		    s->q += 4;
		    memcpy(s->q, s->pkblob, s->pklen);
		    s->q += s->pklen;
		    PUT_32BIT(s->q, s->siglen);
		    s->q += 4;
		    /* Now the data to be signed... */
		    if (!(ssh->remote_bugs & BUG_SSH2_PK_SESSIONID)) {
			PUT_32BIT(s->q, ssh->v2_session_id_len);
			s->q += 4;
		    }
		    memcpy(s->q, ssh->v2_session_id,
			   ssh->v2_session_id_len);
		    s->q += ssh->v2_session_id_len;
		    memcpy(s->q, s->pktout->data + 5,
			   s->pktout->length - 5);
		    s->q += s->pktout->length - 5;
		    /* And finally the (zero) flags word. */
		    PUT_32BIT(s->q, 0);
		    if (!agent_query(s->agentreq, s->len + 4,
				     &vret, &s->retlen,
				     ssh_agent_callback, ssh)) {
			do {
			    crReturnV;
			    if (pktin) {
				bombout(("Unexpected data from server"
					 " while waiting for agent"
					 " response"));
				crStopV;
			    }
			} while (pktin || inlen > 0);
			vret = ssh->agent_response;
			s->retlen = ssh->agent_response_len;
		    }
		    s->ret = vret;
		    sfree(s->agentreq);
		    if (s->ret) {
			if (s->ret[4] == SSH2_AGENT_SIGN_RESPONSE) {
			    logevent("Sending Pageant's response");
			    ssh2_add_sigblob(ssh, s->pktout,
					     s->pkblob, s->pklen,
					     s->ret + 9,
					     GET_32BIT(s->ret + 5));
			    ssh2_pkt_send(ssh, s->pktout);
			    s->type = AUTH_TYPE_PUBLICKEY;
			} else {
			    /* FIXME: less drastic response */
			    bombout(("Pageant failed to answer challenge"));
			    crStopV;
			}
		    }
		}

		/* Do we have any keys left to try? */
		if (s->pkblob_in_agent) {
		    s->done_agent = TRUE;
		    s->tried_pubkey_config = TRUE;
		} else {
		    s->keyi++;
		    if (s->keyi >= s->nkeys)
			s->done_agent = TRUE;
		}

	    } else if (s->can_pubkey && s->publickey_blob &&
		       !s->tried_pubkey_config) {

		struct ssh2_userkey *key;   /* not live over crReturn */
		char *passphrase;	    /* not live over crReturn */

		ssh->pkt_actx = SSH2_PKTCTX_PUBLICKEY;

		s->tried_pubkey_config = TRUE;

		/*
		 * Try the public key supplied in the configuration.
		 *
		 * First, offer the public blob to see if the server is
		 * willing to accept it.
		 */
		s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_REQUEST);
		ssh2_pkt_addstring(s->pktout, s->username);
		ssh2_pkt_addstring(s->pktout, "ssh-connection");
						/* service requested */
		ssh2_pkt_addstring(s->pktout, "publickey");	/* method */
		ssh2_pkt_addbool(s->pktout, FALSE);
						/* no signature included */
		ssh2_pkt_addstring(s->pktout, s->publickey_algorithm);
		ssh2_pkt_addstring_start(s->pktout);
		ssh2_pkt_addstring_data(s->pktout,
					(char *)s->publickey_blob,
					s->publickey_bloblen);
		ssh2_pkt_send(ssh, s->pktout);
		logevent("Offered public key");

		crWaitUntilV(pktin);
		if (pktin->type != SSH2_MSG_USERAUTH_PK_OK) {
		    /* Key refused. Give up. */
		    s->gotit = TRUE; /* reconsider message next loop */
		    s->type = AUTH_TYPE_PUBLICKEY_OFFER_LOUD;
		    continue; /* process this new message */
		}
		logevent("Offer of public key accepted");

		/*
		 * Actually attempt a serious authentication using
		 * the key.
		 */
		if (flags & FLAG_VERBOSE) {
		    c_write_str(ssh, "Authenticating with public key \"");
		    c_write_str(ssh, s->publickey_comment);
		    c_write_str(ssh, "\"\r\n");
		}
		key = NULL;
		while (!key) {
		    const char *error;  /* not live over crReturn */
		    if (s->publickey_encrypted) {
			/*
			 * Get a passphrase from the user.
			 */
			int ret; /* need not be kept over crReturn */
			s->cur_prompt = new_prompts(ssh->frontend);
			s->cur_prompt->to_server = FALSE;
			s->cur_prompt->name = dupstr("SSH key passphrase");
			add_prompt(s->cur_prompt,
				   dupprintf("Passphrase for key \"%.100s\": ",
					     s->publickey_comment),
				   FALSE, SSH_MAX_PASSWORD_LEN);
			ret = get_userpass_input(s->cur_prompt, NULL, 0);
			while (ret < 0) {
			    ssh->send_ok = 1;
			    crWaitUntilV(!pktin);
			    ret = get_userpass_input(s->cur_prompt,
						     in, inlen);
			    ssh->send_ok = 0;
			}
			if (!ret) {
			    /* Failed to get a passphrase. Terminate. */
			    free_prompts(s->cur_prompt);
			    ssh_disconnect(ssh, NULL,
					   "Unable to authenticate",
					   SSH2_DISCONNECT_AUTH_CANCELLED_BY_USER,
					   TRUE);
			    crStopV;
			}
			passphrase =
			    dupstr(s->cur_prompt->prompts[0]->result);
			free_prompts(s->cur_prompt);
		    } else {
			passphrase = NULL; /* no passphrase needed */
		    }

		    /*
		     * Try decrypting the key.
		     */
		    key = ssh2_load_userkey(&ssh->cfg.keyfile, passphrase,
					    &error);
		    if (passphrase) {
			/* burn the evidence */
			memset(passphrase, 0, strlen(passphrase));
			sfree(passphrase);
		    }
		    if (key == SSH2_WRONG_PASSPHRASE || key == NULL) {
			if (passphrase &&
			    (key == SSH2_WRONG_PASSPHRASE)) {
			    c_write_str(ssh, "Wrong passphrase\r\n");
			    key = NULL;
			    /* and loop again */
			} else {
			    c_write_str(ssh, "Unable to load private key (");
			    c_write_str(ssh, error);
			    c_write_str(ssh, ")\r\n");
			    key = NULL;
			    break; /* try something else */
			}
		    }
		}

		if (key) {
		    unsigned char *pkblob, *sigblob, *sigdata;
		    int pkblob_len, sigblob_len, sigdata_len;
		    int p;

		    /*
		     * We have loaded the private key and the server
		     * has announced that it's willing to accept it.
		     * Hallelujah. Generate a signature and send it.
		     */
		    s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_REQUEST);
		    ssh2_pkt_addstring(s->pktout, s->username);
		    ssh2_pkt_addstring(s->pktout, "ssh-connection");
						    /* service requested */
		    ssh2_pkt_addstring(s->pktout, "publickey");
						    /* method */
		    ssh2_pkt_addbool(s->pktout, TRUE);
						    /* signature follows */
		    ssh2_pkt_addstring(s->pktout, key->alg->name);
		    pkblob = key->alg->public_blob(key->data,
						   &pkblob_len);
		    ssh2_pkt_addstring_start(s->pktout);
		    ssh2_pkt_addstring_data(s->pktout, (char *)pkblob,
					    pkblob_len);

		    /*
		     * The data to be signed is:
		     *
		     *   string  session-id
		     *
		     * followed by everything so far placed in the
		     * outgoing packet.
		     */
		    sigdata_len = s->pktout->length - 5 + 4 +
			ssh->v2_session_id_len;
		    if (ssh->remote_bugs & BUG_SSH2_PK_SESSIONID)
			sigdata_len -= 4;
		    sigdata = snewn(sigdata_len, unsigned char);
		    p = 0;
		    if (!(ssh->remote_bugs & BUG_SSH2_PK_SESSIONID)) {
			PUT_32BIT(sigdata+p, ssh->v2_session_id_len);
			p += 4;
		    }
		    memcpy(sigdata+p, ssh->v2_session_id,
			   ssh->v2_session_id_len);
		    p += ssh->v2_session_id_len;
		    memcpy(sigdata+p, s->pktout->data + 5,
			   s->pktout->length - 5);
		    p += s->pktout->length - 5;
		    assert(p == sigdata_len);
		    sigblob = key->alg->sign(key->data, (char *)sigdata,
					     sigdata_len, &sigblob_len);
		    ssh2_add_sigblob(ssh, s->pktout, pkblob, pkblob_len,
				     sigblob, sigblob_len);
		    sfree(pkblob);
		    sfree(sigblob);
		    sfree(sigdata);

		    ssh2_pkt_send(ssh, s->pktout);
                    logevent("Sent public key signature");
		    s->type = AUTH_TYPE_PUBLICKEY;
		    key->alg->freekey(key->data);
		}

#ifndef NO_GSSAPI
	    } else if (s->can_gssapi && !s->tried_gssapi) {

		/* GSSAPI Authentication */

		int micoffset, len;
		char *data;
		Ssh_gss_buf mic;
		s->type = AUTH_TYPE_GSSAPI;
		s->tried_gssapi = TRUE;
		s->gotit = TRUE;
		ssh->pkt_actx = SSH2_PKTCTX_GSSAPI;

		/*
		 * Pick the highest GSS library on the preference
		 * list.
		 */
		{
		    int i, j;
		    s->gsslib = NULL;
		    for (i = 0; i < ngsslibs; i++) {
			int want_id = ssh->cfg.ssh_gsslist[i];
			for (j = 0; j < ssh->gsslibs->nlibraries; j++)
			    if (ssh->gsslibs->libraries[j].id == want_id) {
				s->gsslib = &ssh->gsslibs->libraries[j];
				goto got_gsslib;   /* double break */
			    }
		    }
		    got_gsslib:
		    /*
		     * We always expect to have found something in
		     * the above loop: we only came here if there
		     * was at least one viable GSS library, and the
		     * preference list should always mention
		     * everything and only change the order.
		     */
		    assert(s->gsslib);
		}

		if (s->gsslib->gsslogmsg)
		    logevent(s->gsslib->gsslogmsg);

		/* Sending USERAUTH_REQUEST with "gssapi-with-mic" method */
		s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_REQUEST);
		ssh2_pkt_addstring(s->pktout, s->username);
		ssh2_pkt_addstring(s->pktout, "ssh-connection");
		ssh2_pkt_addstring(s->pktout, "gssapi-with-mic");
                logevent("Attempting GSSAPI authentication");

		/* add mechanism info */
		s->gsslib->indicate_mech(s->gsslib, &s->gss_buf);

		/* number of GSSAPI mechanisms */
		ssh2_pkt_adduint32(s->pktout,1);

		/* length of OID + 2 */
		ssh2_pkt_adduint32(s->pktout, s->gss_buf.length + 2);
		ssh2_pkt_addbyte(s->pktout, SSH2_GSS_OIDTYPE);

		/* length of OID */
		ssh2_pkt_addbyte(s->pktout, (unsigned char) s->gss_buf.length);

		ssh_pkt_adddata(s->pktout, s->gss_buf.value,
				s->gss_buf.length);
		ssh2_pkt_send(ssh, s->pktout);
		crWaitUntilV(pktin);
		if (pktin->type != SSH2_MSG_USERAUTH_GSSAPI_RESPONSE) {
		    logevent("GSSAPI authentication request refused");
		    continue;
		}

		/* check returned packet ... */

		ssh_pkt_getstring(pktin, &data, &len);
		s->gss_rcvtok.value = data;
		s->gss_rcvtok.length = len;
		if (s->gss_rcvtok.length != s->gss_buf.length + 2 ||
		    ((char *)s->gss_rcvtok.value)[0] != SSH2_GSS_OIDTYPE ||
		    ((char *)s->gss_rcvtok.value)[1] != s->gss_buf.length ||
		    memcmp((char *)s->gss_rcvtok.value + 2,
			   s->gss_buf.value,s->gss_buf.length) ) {
		    logevent("GSSAPI authentication - wrong response from server");
		    continue;
		}

		/* now start running */
		s->gss_stat = s->gsslib->import_name(s->gsslib,
						     ssh->fullhostname,
						     &s->gss_srv_name);
		if (s->gss_stat != SSH_GSS_OK) {
		    if (s->gss_stat == SSH_GSS_BAD_HOST_NAME)
			logevent("GSSAPI import name failed - Bad service name");
		    else
			logevent("GSSAPI import name failed");
		    continue;
		}

		/* fetch TGT into GSS engine */
		s->gss_stat = s->gsslib->acquire_cred(s->gsslib, &s->gss_ctx);

		if (s->gss_stat != SSH_GSS_OK) {
		    logevent("GSSAPI authentication failed to get credentials");
		    s->gsslib->release_name(s->gsslib, &s->gss_srv_name);
		    continue;
		}

		/* initial tokens are empty */
		SSH_GSS_CLEAR_BUF(&s->gss_rcvtok);
		SSH_GSS_CLEAR_BUF(&s->gss_sndtok);

		/* now enter the loop */
		do {
		    s->gss_stat = s->gsslib->init_sec_context
			(s->gsslib,
			 &s->gss_ctx,
			 s->gss_srv_name,
			 ssh->cfg.gssapifwd,
			 &s->gss_rcvtok,
			 &s->gss_sndtok);

		    if (s->gss_stat!=SSH_GSS_S_COMPLETE &&
			s->gss_stat!=SSH_GSS_S_CONTINUE_NEEDED) {
			logevent("GSSAPI authentication initialisation failed");

			if (s->gsslib->display_status(s->gsslib, s->gss_ctx,
						      &s->gss_buf) == SSH_GSS_OK) {
			    logevent(s->gss_buf.value);
			    sfree(s->gss_buf.value);
			}

			break;
		    }
		    logevent("GSSAPI authentication initialised");

		    /* Client and server now exchange tokens until GSSAPI
		     * no longer says CONTINUE_NEEDED */

		    if (s->gss_sndtok.length != 0) {
			s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_GSSAPI_TOKEN);
			ssh_pkt_addstring_start(s->pktout);
			ssh_pkt_addstring_data(s->pktout,s->gss_sndtok.value,s->gss_sndtok.length);
			ssh2_pkt_send(ssh, s->pktout);
			s->gsslib->free_tok(s->gsslib, &s->gss_sndtok);
		    }

		    if (s->gss_stat == SSH_GSS_S_CONTINUE_NEEDED) {
			crWaitUntilV(pktin);
			if (pktin->type != SSH2_MSG_USERAUTH_GSSAPI_TOKEN) {
			    logevent("GSSAPI authentication - bad server response");
			    s->gss_stat = SSH_GSS_FAILURE;
			    break;
			}
			ssh_pkt_getstring(pktin, &data, &len);
			s->gss_rcvtok.value = data;
			s->gss_rcvtok.length = len;
		    }
		} while (s-> gss_stat == SSH_GSS_S_CONTINUE_NEEDED);

		if (s->gss_stat != SSH_GSS_OK) {
		    s->gsslib->release_name(s->gsslib, &s->gss_srv_name);
		    s->gsslib->release_cred(s->gsslib, &s->gss_ctx);
		    continue;
		}
		logevent("GSSAPI authentication loop finished OK");

		/* Now send the MIC */

		s->pktout = ssh2_pkt_init(0);
		micoffset = s->pktout->length;
		ssh_pkt_addstring_start(s->pktout);
		ssh_pkt_addstring_data(s->pktout, (char *)ssh->v2_session_id, ssh->v2_session_id_len);
		ssh_pkt_addbyte(s->pktout, SSH2_MSG_USERAUTH_REQUEST);
		ssh_pkt_addstring(s->pktout, s->username);
		ssh_pkt_addstring(s->pktout, "ssh-connection");
		ssh_pkt_addstring(s->pktout, "gssapi-with-mic");

		s->gss_buf.value = (char *)s->pktout->data + micoffset;
		s->gss_buf.length = s->pktout->length - micoffset;

		s->gsslib->get_mic(s->gsslib, s->gss_ctx, &s->gss_buf, &mic);
		s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_GSSAPI_MIC);
		ssh_pkt_addstring_start(s->pktout);
		ssh_pkt_addstring_data(s->pktout, mic.value, mic.length);
		ssh2_pkt_send(ssh, s->pktout);
		s->gsslib->free_mic(s->gsslib, &mic);

		s->gotit = FALSE;

		s->gsslib->release_name(s->gsslib, &s->gss_srv_name);
		s->gsslib->release_cred(s->gsslib, &s->gss_ctx);
		continue;
#endif
	    } else if (s->can_keyb_inter && !s->kbd_inter_refused) {

		/*
		 * Keyboard-interactive authentication.
		 */

		s->type = AUTH_TYPE_KEYBOARD_INTERACTIVE;

		ssh->pkt_actx = SSH2_PKTCTX_KBDINTER;

		s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_REQUEST);
		ssh2_pkt_addstring(s->pktout, s->username);
		ssh2_pkt_addstring(s->pktout, "ssh-connection");
							/* service requested */
		ssh2_pkt_addstring(s->pktout, "keyboard-interactive");
							/* method */
		ssh2_pkt_addstring(s->pktout, "");	/* lang */
		ssh2_pkt_addstring(s->pktout, "");	/* submethods */
		ssh2_pkt_send(ssh, s->pktout);
                
                logevent("Attempting keyboard-interactive authentication");

		crWaitUntilV(pktin);
		if (pktin->type != SSH2_MSG_USERAUTH_INFO_REQUEST) {
		    /* Server is not willing to do keyboard-interactive
		     * at all (or, bizarrely but legally, accepts the
		     * user without actually issuing any prompts).
		     * Give up on it entirely. */
		    s->gotit = TRUE;
		    s->type = AUTH_TYPE_KEYBOARD_INTERACTIVE_QUIET;
		    s->kbd_inter_refused = TRUE; /* don't try it again */
		    continue;
		}

		/*
		 * Loop while the server continues to send INFO_REQUESTs.
		 */
		while (pktin->type == SSH2_MSG_USERAUTH_INFO_REQUEST) {

		    char *name, *inst, *lang;
		    int name_len, inst_len, lang_len;
		    int i;

		    /*
		     * We've got a fresh USERAUTH_INFO_REQUEST.
		     * Get the preamble and start building a prompt.
		     */
		    ssh_pkt_getstring(pktin, &name, &name_len);
		    ssh_pkt_getstring(pktin, &inst, &inst_len);
		    ssh_pkt_getstring(pktin, &lang, &lang_len);
		    s->cur_prompt = new_prompts(ssh->frontend);
		    s->cur_prompt->to_server = TRUE;

		    /*
		     * Get any prompt(s) from the packet.
		     */
		    s->num_prompts = ssh_pkt_getuint32(pktin);
		    for (i = 0; i < s->num_prompts; i++) {
			char *prompt;
			int prompt_len;
			int echo;
			static char noprompt[] =
			    "<server failed to send prompt>: ";

			ssh_pkt_getstring(pktin, &prompt, &prompt_len);
			echo = ssh2_pkt_getbool(pktin);
			if (!prompt_len) {
			    prompt = noprompt;
			    prompt_len = lenof(noprompt)-1;
			}
			add_prompt(s->cur_prompt,
				   dupprintf("%.*s", prompt_len, prompt),
				   echo, SSH_MAX_PASSWORD_LEN);
		    }

		    if (name_len) {
			/* FIXME: better prefix to distinguish from
			 * local prompts? */
			s->cur_prompt->name =
			    dupprintf("SSH server: %.*s", name_len, name);
			s->cur_prompt->name_reqd = TRUE;
		    } else {
			s->cur_prompt->name =
			    dupstr("SSH server authentication");
			s->cur_prompt->name_reqd = FALSE;
		    }
		    /* We add a prefix to try to make it clear that a prompt
		     * has come from the server.
		     * FIXME: ugly to print "Using..." in prompt _every_
		     * time round. Can this be done more subtly? */
		    /* Special case: for reasons best known to themselves,
		     * some servers send k-i requests with no prompts and
		     * nothing to display. Keep quiet in this case. */
		    if (s->num_prompts || name_len || inst_len) {
			s->cur_prompt->instruction =
			    dupprintf("Using keyboard-interactive authentication.%s%.*s",
				      inst_len ? "\n" : "", inst_len, inst);
			s->cur_prompt->instr_reqd = TRUE;
		    } else {
			s->cur_prompt->instr_reqd = FALSE;
		    }

		    /*
                     * Display any instructions, and get the user's
                     * response(s).
		     */
		    {
			int ret; /* not live over crReturn */
			ret = get_userpass_input(s->cur_prompt, NULL, 0);
			while (ret < 0) {
			    ssh->send_ok = 1;
			    crWaitUntilV(!pktin);
			    ret = get_userpass_input(s->cur_prompt, in, inlen);
			    ssh->send_ok = 0;
			}
			if (!ret) {
			    /*
			     * Failed to get responses. Terminate.
			     */
			    free_prompts(s->cur_prompt);
			    ssh_disconnect(ssh, NULL, "Unable to authenticate",
					   SSH2_DISCONNECT_AUTH_CANCELLED_BY_USER,
					   TRUE);
			    crStopV;
			}
		    }

		    /*
		     * Send the response(s) to the server.
		     */
		    s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_INFO_RESPONSE);
		    ssh2_pkt_adduint32(s->pktout, s->num_prompts);
		    for (i=0; i < s->num_prompts; i++) {
			dont_log_password(ssh, s->pktout, PKTLOG_BLANK);
			ssh2_pkt_addstring(s->pktout,
					   s->cur_prompt->prompts[i]->result);
			end_log_omission(ssh, s->pktout);
		    }
		    ssh2_pkt_send_with_padding(ssh, s->pktout, 256);

                    /*
                     * Free the prompts structure from this iteration.
                     * If there's another, a new one will be allocated
                     * when we return to the top of this while loop.
                     */
                    free_prompts(s->cur_prompt);

		    /*
		     * Get the next packet in case it's another
		     * INFO_REQUEST.
		     */
		    crWaitUntilV(pktin);

		}

		/*
		 * We should have SUCCESS or FAILURE now.
		 */
		s->gotit = TRUE;

	    } else if (s->can_passwd) {

		/*
		 * Plain old password authentication.
		 */
		int ret; /* not live over crReturn */
		int changereq_first_time; /* not live over crReturn */

		ssh->pkt_actx = SSH2_PKTCTX_PASSWORD;

		s->cur_prompt = new_prompts(ssh->frontend);
		s->cur_prompt->to_server = TRUE;
		s->cur_prompt->name = dupstr("SSH password");
		add_prompt(s->cur_prompt, dupprintf("%.90s@%.90s's password: ",
						    s->username,
						    ssh->savedhost),
			   FALSE, SSH_MAX_PASSWORD_LEN);

		ret = get_userpass_input(s->cur_prompt, NULL, 0);
		while (ret < 0) {
		    ssh->send_ok = 1;
		    crWaitUntilV(!pktin);
		    ret = get_userpass_input(s->cur_prompt, in, inlen);
		    ssh->send_ok = 0;
		}
		if (!ret) {
		    /*
		     * Failed to get responses. Terminate.
		     */
		    free_prompts(s->cur_prompt);
		    ssh_disconnect(ssh, NULL, "Unable to authenticate",
				   SSH2_DISCONNECT_AUTH_CANCELLED_BY_USER,
				   TRUE);
		    crStopV;
		}
		/*
		 * Squirrel away the password. (We may need it later if
		 * asked to change it.)
		 */
		s->password = dupstr(s->cur_prompt->prompts[0]->result);
		free_prompts(s->cur_prompt);

		/*
		 * Send the password packet.
		 *
		 * We pad out the password packet to 256 bytes to make
		 * it harder for an attacker to find the length of the
		 * user's password.
		 *
		 * Anyone using a password longer than 256 bytes
		 * probably doesn't have much to worry about from
		 * people who find out how long their password is!
		 */
		s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_REQUEST);
		ssh2_pkt_addstring(s->pktout, s->username);
		ssh2_pkt_addstring(s->pktout, "ssh-connection");
							/* service requested */
		ssh2_pkt_addstring(s->pktout, "password");
		ssh2_pkt_addbool(s->pktout, FALSE);
		dont_log_password(ssh, s->pktout, PKTLOG_BLANK);
		ssh2_pkt_addstring(s->pktout, s->password);
		end_log_omission(ssh, s->pktout);
		ssh2_pkt_send_with_padding(ssh, s->pktout, 256);
		logevent("Sent password");
		s->type = AUTH_TYPE_PASSWORD;

		/*
		 * Wait for next packet, in case it's a password change
		 * request.
		 */
		crWaitUntilV(pktin);
		changereq_first_time = TRUE;

		while (pktin->type == SSH2_MSG_USERAUTH_PASSWD_CHANGEREQ) {

		    /* 
		     * We're being asked for a new password
		     * (perhaps not for the first time).
		     * Loop until the server accepts it.
		     */

		    int got_new = FALSE; /* not live over crReturn */
		    char *prompt;   /* not live over crReturn */
		    int prompt_len; /* not live over crReturn */
		    
		    {
			char *msg;
			if (changereq_first_time)
			    msg = "Server requested password change";
			else
			    msg = "Server rejected new password";
			logevent(msg);
			c_write_str(ssh, msg);
			c_write_str(ssh, "\r\n");
		    }

		    ssh_pkt_getstring(pktin, &prompt, &prompt_len);

		    s->cur_prompt = new_prompts(ssh->frontend);
		    s->cur_prompt->to_server = TRUE;
		    s->cur_prompt->name = dupstr("New SSH password");
		    s->cur_prompt->instruction =
			dupprintf("%.*s", prompt_len, prompt);
		    s->cur_prompt->instr_reqd = TRUE;
		    /*
		     * There's no explicit requirement in the protocol
		     * for the "old" passwords in the original and
		     * password-change messages to be the same, and
		     * apparently some Cisco kit supports password change
		     * by the user entering a blank password originally
		     * and the real password subsequently, so,
		     * reluctantly, we prompt for the old password again.
		     *
		     * (On the other hand, some servers don't even bother
		     * to check this field.)
		     */
		    add_prompt(s->cur_prompt,
			       dupstr("Current password (blank for previously entered password): "),
			       FALSE, SSH_MAX_PASSWORD_LEN);
		    add_prompt(s->cur_prompt, dupstr("Enter new password: "),
			       FALSE, SSH_MAX_PASSWORD_LEN);
		    add_prompt(s->cur_prompt, dupstr("Confirm new password: "),
			       FALSE, SSH_MAX_PASSWORD_LEN);

		    /*
		     * Loop until the user manages to enter the same
		     * password twice.
		     */
		    while (!got_new) {

			ret = get_userpass_input(s->cur_prompt, NULL, 0);
			while (ret < 0) {
			    ssh->send_ok = 1;
			    crWaitUntilV(!pktin);
			    ret = get_userpass_input(s->cur_prompt, in, inlen);
			    ssh->send_ok = 0;
			}
			if (!ret) {
			    /*
			     * Failed to get responses. Terminate.
			     */
			    /* burn the evidence */
			    free_prompts(s->cur_prompt);
			    memset(s->password, 0, strlen(s->password));
			    sfree(s->password);
			    ssh_disconnect(ssh, NULL, "Unable to authenticate",
					   SSH2_DISCONNECT_AUTH_CANCELLED_BY_USER,
					   TRUE);
			    crStopV;
			}

			/*
			 * If the user specified a new original password
			 * (IYSWIM), overwrite any previously specified
			 * one.
			 * (A side effect is that the user doesn't have to
			 * re-enter it if they louse up the new password.)
			 */
			if (s->cur_prompt->prompts[0]->result[0]) {
			    memset(s->password, 0, strlen(s->password));
				/* burn the evidence */
			    sfree(s->password);
			    s->password =
				dupstr(s->cur_prompt->prompts[0]->result);
			}

			/*
			 * Check the two new passwords match.
			 */
			got_new = (strcmp(s->cur_prompt->prompts[1]->result,
					  s->cur_prompt->prompts[2]->result)
				   == 0);
			if (!got_new)
			    /* They don't. Silly user. */
			    c_write_str(ssh, "Passwords do not match\r\n");

		    }

		    /*
		     * Send the new password (along with the old one).
		     * (see above for padding rationale)
		     */
		    s->pktout = ssh2_pkt_init(SSH2_MSG_USERAUTH_REQUEST);
		    ssh2_pkt_addstring(s->pktout, s->username);
		    ssh2_pkt_addstring(s->pktout, "ssh-connection");
							/* service requested */
		    ssh2_pkt_addstring(s->pktout, "password");
		    ssh2_pkt_addbool(s->pktout, TRUE);
		    dont_log_password(ssh, s->pktout, PKTLOG_BLANK);
		    ssh2_pkt_addstring(s->pktout, s->password);
		    ssh2_pkt_addstring(s->pktout,
				       s->cur_prompt->prompts[1]->result);
		    free_prompts(s->cur_prompt);
		    end_log_omission(ssh, s->pktout);
		    ssh2_pkt_send_with_padding(ssh, s->pktout, 256);
		    logevent("Sent new password");
		    
		    /*
		     * Now see what the server has to say about it.
		     * (If it's CHANGEREQ again, it's not happy with the
		     * new password.)
		     */
		    crWaitUntilV(pktin);
		    changereq_first_time = FALSE;

		}

		/*
		 * We need to reexamine the current pktin at the top
		 * of the loop. Either:
		 *  - we weren't asked to change password at all, in
		 *    which case it's a SUCCESS or FAILURE with the
		 *    usual meaning
		 *  - we sent a new password, and the server was
		 *    either OK with it (SUCCESS or FAILURE w/partial
		 *    success) or unhappy with the _old_ password
		 *    (FAILURE w/o partial success)
		 * In any of these cases, we go back to the top of
		 * the loop and start again.
		 */
		s->gotit = TRUE;

		/*
		 * We don't need the old password any more, in any
		 * case. Burn the evidence.
		 */
		memset(s->password, 0, strlen(s->password));
		sfree(s->password);

	    } else {
		char *str = dupprintf("No supported authentication methods available"
				      " (server sent: %.*s)",
				      methlen, methods);

		ssh_disconnect(ssh, str,
			       "No supported authentication methods available",
			       SSH2_DISCONNECT_NO_MORE_AUTH_METHODS_AVAILABLE,
			       FALSE);
		sfree(str);

		crStopV;

	    }

	}
    }
    ssh->packet_dispatch[SSH2_MSG_USERAUTH_BANNER] = NULL;

    /* Clear up various bits and pieces from authentication. */
    if (s->publickey_blob) {
	sfree(s->publickey_blob);
	sfree(s->publickey_comment);
    }
    if (s->agent_response)
	sfree(s->agent_response);

    if (s->userauth_success) {
	/*
	 * We've just received USERAUTH_SUCCESS, and we haven't sent any
	 * packets since. Signal the transport layer to consider enacting
	 * delayed compression.
	 *
	 * (Relying on we_are_in is not sufficient, as
	 * draft-miller-secsh-compression-delayed is quite clear that it
	 * triggers on USERAUTH_SUCCESS specifically, and we_are_in can
	 * become set for other reasons.)
	 */
	do_ssh2_transport(ssh, "enabling delayed compression", -2, NULL);
    }

    /*
     * Now the connection protocol has started, one way or another.
     */

    ssh->channels = newtree234(ssh_channelcmp);

    /*
     * Set up handlers for some connection protocol messages, so we
     * don't have to handle them repeatedly in this coroutine.
     */
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_WINDOW_ADJUST] =
	ssh2_msg_channel_window_adjust;
    ssh->packet_dispatch[SSH2_MSG_GLOBAL_REQUEST] =
	ssh2_msg_global_request;

    /*
     * Create the main session channel.
     */
    if (ssh->cfg.ssh_no_shell) {
	ssh->mainchan = NULL;
    } else if (*ssh->cfg.ssh_nc_host) {
	/*
	 * Just start a direct-tcpip channel and use it as the main
	 * channel.
	 */
	ssh->mainchan = snew(struct ssh_channel);
	ssh->mainchan->ssh = ssh;
	ssh2_channel_init(ssh->mainchan);
	logeventf(ssh,
		  "Opening direct-tcpip channel to %s:%d in place of session",
		  ssh->cfg.ssh_nc_host, ssh->cfg.ssh_nc_port);
	s->pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_OPEN);
	ssh2_pkt_addstring(s->pktout, "direct-tcpip");
	ssh2_pkt_adduint32(s->pktout, ssh->mainchan->localid);
	ssh2_pkt_adduint32(s->pktout, ssh->mainchan->v.v2.locwindow);/* our window size */
	ssh2_pkt_adduint32(s->pktout, OUR_V2_MAXPKT);      /* our max pkt size */
	ssh2_pkt_addstring(s->pktout, ssh->cfg.ssh_nc_host);
	ssh2_pkt_adduint32(s->pktout, ssh->cfg.ssh_nc_port);
	/*
	 * There's nothing meaningful to put in the originator
	 * fields, but some servers insist on syntactically correct
	 * information.
	 */
	ssh2_pkt_addstring(s->pktout, "0.0.0.0");
	ssh2_pkt_adduint32(s->pktout, 0);
	ssh2_pkt_send(ssh, s->pktout);

	crWaitUntilV(pktin);
	if (pktin->type != SSH2_MSG_CHANNEL_OPEN_CONFIRMATION) {
	    bombout(("Server refused to open a direct-tcpip channel"));
	    crStopV;
	    /* FIXME: error data comes back in FAILURE packet */
	}
	if (ssh_pkt_getuint32(pktin) != ssh->mainchan->localid) {
	    bombout(("Server's channel confirmation cited wrong channel"));
	    crStopV;
	}
	ssh->mainchan->remoteid = ssh_pkt_getuint32(pktin);
	ssh->mainchan->halfopen = FALSE;
	ssh->mainchan->type = CHAN_MAINSESSION;
	ssh->mainchan->v.v2.remwindow = ssh_pkt_getuint32(pktin);
	ssh->mainchan->v.v2.remmaxpkt = ssh_pkt_getuint32(pktin);
	add234(ssh->channels, ssh->mainchan);
	update_specials_menu(ssh->frontend);
	logevent("Opened direct-tcpip channel");
	ssh->ncmode = TRUE;
    } else {
	ssh->mainchan = snew(struct ssh_channel);
	ssh->mainchan->ssh = ssh;
	ssh2_channel_init(ssh->mainchan);
	s->pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_OPEN);
	ssh2_pkt_addstring(s->pktout, "session");
	ssh2_pkt_adduint32(s->pktout, ssh->mainchan->localid);
	ssh2_pkt_adduint32(s->pktout, ssh->mainchan->v.v2.locwindow);/* our window size */
	ssh2_pkt_adduint32(s->pktout, OUR_V2_MAXPKT);    /* our max pkt size */
	ssh2_pkt_send(ssh, s->pktout);
	crWaitUntilV(pktin);
	if (pktin->type != SSH2_MSG_CHANNEL_OPEN_CONFIRMATION) {
	    bombout(("Server refused to open a session"));
	    crStopV;
	    /* FIXME: error data comes back in FAILURE packet */
	}
	if (ssh_pkt_getuint32(pktin) != ssh->mainchan->localid) {
	    bombout(("Server's channel confirmation cited wrong channel"));
	    crStopV;
	}
	ssh->mainchan->remoteid = ssh_pkt_getuint32(pktin);
	ssh->mainchan->halfopen = FALSE;
	ssh->mainchan->type = CHAN_MAINSESSION;
	ssh->mainchan->v.v2.remwindow = ssh_pkt_getuint32(pktin);
	ssh->mainchan->v.v2.remmaxpkt = ssh_pkt_getuint32(pktin);
	add234(ssh->channels, ssh->mainchan);
	update_specials_menu(ssh->frontend);
	logevent("Opened channel for session");
	ssh->ncmode = FALSE;
    }

    /*
     * Now we have a channel, make dispatch table entries for
     * general channel-based messages.
     */
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_DATA] =
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_EXTENDED_DATA] =
	ssh2_msg_channel_data;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_EOF] = ssh2_msg_channel_eof;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_CLOSE] = ssh2_msg_channel_close;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_OPEN_CONFIRMATION] =
	ssh2_msg_channel_open_confirmation;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_OPEN_FAILURE] =
	ssh2_msg_channel_open_failure;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_REQUEST] =
	ssh2_msg_channel_request;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_OPEN] =
	ssh2_msg_channel_open;

    if (ssh->mainchan && ssh->cfg.ssh_simple) {
	/*
	 * This message indicates to the server that we promise
	 * not to try to run any other channel in parallel with
	 * this one, so it's safe for it to advertise a very large
	 * window and leave the flow control to TCP.
	 */
	s->pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
	ssh2_pkt_adduint32(s->pktout, ssh->mainchan->remoteid);
	ssh2_pkt_addstring(s->pktout, "simple@putty.projects.tartarus.org");
	ssh2_pkt_addbool(s->pktout, 0); /* no reply */
	ssh2_pkt_send(ssh, s->pktout);
    }

    /*
     * Potentially enable X11 forwarding.
     */
    if (ssh->mainchan && !ssh->ncmode && ssh->cfg.x11_forward &&
	(ssh->x11disp = x11_setup_display(ssh->cfg.x11_display,
					  ssh->cfg.x11_auth, &ssh->cfg))) {
	logevent("Requesting X11 forwarding");
	s->pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
	ssh2_pkt_adduint32(s->pktout, ssh->mainchan->remoteid);
	ssh2_pkt_addstring(s->pktout, "x11-req");
	ssh2_pkt_addbool(s->pktout, 1);	       /* want reply */
	ssh2_pkt_addbool(s->pktout, 0);	       /* many connections */
	ssh2_pkt_addstring(s->pktout, ssh->x11disp->remoteauthprotoname);
	/*
	 * Note that while we blank the X authentication data here, we don't
	 * take any special action to blank the start of an X11 channel,
	 * so using MIT-MAGIC-COOKIE-1 and actually opening an X connection
	 * without having session blanking enabled is likely to leak your
	 * cookie into the log.
	 */
	dont_log_password(ssh, s->pktout, PKTLOG_BLANK);
	ssh2_pkt_addstring(s->pktout, ssh->x11disp->remoteauthdatastring);
	end_log_omission(ssh, s->pktout);
	ssh2_pkt_adduint32(s->pktout, ssh->x11disp->screennum);
	ssh2_pkt_send(ssh, s->pktout);

	crWaitUntilV(pktin);

	if (pktin->type != SSH2_MSG_CHANNEL_SUCCESS) {
	    if (pktin->type != SSH2_MSG_CHANNEL_FAILURE) {
		bombout(("Unexpected response to X11 forwarding request:"
			 " packet type %d", pktin->type));
		crStopV;
	    }
	    logevent("X11 forwarding refused");
	} else {
	    logevent("X11 forwarding enabled");
	    ssh->X11_fwd_enabled = TRUE;
	}
    }

    /*
     * Enable port forwardings.
     */
    ssh_setup_portfwd(ssh, &ssh->cfg);

    /*
     * Potentially enable agent forwarding.
     */
    if (ssh->mainchan && !ssh->ncmode && ssh->cfg.agentfwd && agent_exists()) {
	logevent("Requesting OpenSSH-style agent forwarding");
	s->pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
	ssh2_pkt_adduint32(s->pktout, ssh->mainchan->remoteid);
	ssh2_pkt_addstring(s->pktout, "auth-agent-req@openssh.com");
	ssh2_pkt_addbool(s->pktout, 1);	       /* want reply */
	ssh2_pkt_send(ssh, s->pktout);

	crWaitUntilV(pktin);

	if (pktin->type != SSH2_MSG_CHANNEL_SUCCESS) {
	    if (pktin->type != SSH2_MSG_CHANNEL_FAILURE) {
		bombout(("Unexpected response to agent forwarding request:"
			 " packet type %d", pktin->type));
		crStopV;
	    }
	    logevent("Agent forwarding refused");
	} else {
	    logevent("Agent forwarding enabled");
	    ssh->agentfwd_enabled = TRUE;
	}
    }

    /*
     * Now allocate a pty for the session.
     */
    if (ssh->mainchan && !ssh->ncmode && !ssh->cfg.nopty) {
	/* Unpick the terminal-speed string. */
	/* XXX perhaps we should allow no speeds to be sent. */
        ssh->ospeed = 38400; ssh->ispeed = 38400; /* last-resort defaults */
	sscanf(ssh->cfg.termspeed, "%d,%d", &ssh->ospeed, &ssh->ispeed);
	/* Build the pty request. */
	s->pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
	ssh2_pkt_adduint32(s->pktout, ssh->mainchan->remoteid);	/* recipient channel */
	ssh2_pkt_addstring(s->pktout, "pty-req");
	ssh2_pkt_addbool(s->pktout, 1);	       /* want reply */
	ssh2_pkt_addstring(s->pktout, ssh->cfg.termtype);
	ssh2_pkt_adduint32(s->pktout, ssh->term_width);
	ssh2_pkt_adduint32(s->pktout, ssh->term_height);
	ssh2_pkt_adduint32(s->pktout, 0);	       /* pixel width */
	ssh2_pkt_adduint32(s->pktout, 0);	       /* pixel height */
	ssh2_pkt_addstring_start(s->pktout);
	parse_ttymodes(ssh, ssh->cfg.ttymodes,
		       ssh2_send_ttymode, (void *)s->pktout);
	ssh2_pkt_addbyte(s->pktout, SSH2_TTY_OP_ISPEED);
	ssh2_pkt_adduint32(s->pktout, ssh->ispeed);
	ssh2_pkt_addbyte(s->pktout, SSH2_TTY_OP_OSPEED);
	ssh2_pkt_adduint32(s->pktout, ssh->ospeed);
	ssh2_pkt_addstring_data(s->pktout, "\0", 1); /* TTY_OP_END */
	ssh2_pkt_send(ssh, s->pktout);
	ssh->state = SSH_STATE_INTERMED;

	crWaitUntilV(pktin);

	if (pktin->type != SSH2_MSG_CHANNEL_SUCCESS) {
	    if (pktin->type != SSH2_MSG_CHANNEL_FAILURE) {
		bombout(("Unexpected response to pty request:"
			 " packet type %d", pktin->type));
		crStopV;
	    }
	    c_write_str(ssh, "Server refused to allocate pty\r\n");
	    ssh->editing = ssh->echoing = 1;
	} else {
	    logeventf(ssh, "Allocated pty (ospeed %dbps, ispeed %dbps)",
		      ssh->ospeed, ssh->ispeed);
	}
    } else {
	ssh->editing = ssh->echoing = 1;
    }

    /*
     * Send environment variables.
     * 
     * Simplest thing here is to send all the requests at once, and
     * then wait for a whole bunch of successes or failures.
     */
    if (ssh->mainchan && !ssh->ncmode && *ssh->cfg.environmt) {
	char *e = ssh->cfg.environmt;
	char *var, *varend, *val;

	s->num_env = 0;

	while (*e) {
	    var = e;
	    while (*e && *e != '\t') e++;
	    varend = e;
	    if (*e == '\t') e++;
	    val = e;
	    while (*e) e++;
	    e++;

	    s->pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
	    ssh2_pkt_adduint32(s->pktout, ssh->mainchan->remoteid);
	    ssh2_pkt_addstring(s->pktout, "env");
	    ssh2_pkt_addbool(s->pktout, 1);	       /* want reply */
	    ssh2_pkt_addstring_start(s->pktout);
	    ssh2_pkt_addstring_data(s->pktout, var, varend-var);
	    ssh2_pkt_addstring(s->pktout, val);
	    ssh2_pkt_send(ssh, s->pktout);

	    s->num_env++;
	}

	logeventf(ssh, "Sent %d environment variables", s->num_env);

	s->env_ok = 0;
	s->env_left = s->num_env;

	while (s->env_left > 0) {
	    crWaitUntilV(pktin);

	    if (pktin->type != SSH2_MSG_CHANNEL_SUCCESS) {
		if (pktin->type != SSH2_MSG_CHANNEL_FAILURE) {
		    bombout(("Unexpected response to environment request:"
			     " packet type %d", pktin->type));
		    crStopV;
		}
	    } else {
		s->env_ok++;
	    }

	    s->env_left--;
	}

	if (s->env_ok == s->num_env) {
	    logevent("All environment variables successfully set");
	} else if (s->env_ok == 0) {
	    logevent("All environment variables refused");
	    c_write_str(ssh, "Server refused to set environment variables\r\n");
	} else {
	    logeventf(ssh, "%d environment variables refused",
		      s->num_env - s->env_ok);
	    c_write_str(ssh, "Server refused to set all environment variables\r\n");
	}
    }

    /*
     * Start a shell or a remote command. We may have to attempt
     * this twice if the config data has provided a second choice
     * of command.
     */
    if (ssh->mainchan && !ssh->ncmode) while (1) {
	int subsys;
	char *cmd;

	if (ssh->fallback_cmd) {
	    subsys = ssh->cfg.ssh_subsys2;
	    cmd = ssh->cfg.remote_cmd_ptr2;
	} else {
	    subsys = ssh->cfg.ssh_subsys;
	    cmd = ssh->cfg.remote_cmd_ptr;
	    if (!cmd) cmd = ssh->cfg.remote_cmd;
	}

	s->pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
	ssh2_pkt_adduint32(s->pktout, ssh->mainchan->remoteid);	/* recipient channel */
	if (subsys) {
	    ssh2_pkt_addstring(s->pktout, "subsystem");
	    ssh2_pkt_addbool(s->pktout, 1);	       /* want reply */
	    ssh2_pkt_addstring(s->pktout, cmd);
	} else if (*cmd) {
	    ssh2_pkt_addstring(s->pktout, "exec");
	    ssh2_pkt_addbool(s->pktout, 1);	       /* want reply */
	    ssh2_pkt_addstring(s->pktout, cmd);
	} else {
	    ssh2_pkt_addstring(s->pktout, "shell");
	    ssh2_pkt_addbool(s->pktout, 1);	       /* want reply */
	}
	ssh2_pkt_send(ssh, s->pktout);

	crWaitUntilV(pktin);

	if (pktin->type != SSH2_MSG_CHANNEL_SUCCESS) {
	    if (pktin->type != SSH2_MSG_CHANNEL_FAILURE) {
		bombout(("Unexpected response to shell/command request:"
			 " packet type %d", pktin->type));
		crStopV;
	    }
	    /*
	     * We failed to start the command. If this is the
	     * fallback command, we really are finished; if it's
	     * not, and if the fallback command exists, try falling
	     * back to it before complaining.
	     */
	    if (!ssh->fallback_cmd && ssh->cfg.remote_cmd_ptr2 != NULL) {
		logevent("Primary command failed; attempting fallback");
		ssh->fallback_cmd = TRUE;
		continue;
	    }
	    bombout(("Server refused to start a shell/command"));
	    crStopV;
	} else {
	    logevent("Started a shell/command");
	}
	break;
    }

    ssh->state = SSH_STATE_SESSION;
    if (ssh->size_needed)
	ssh_size(ssh, ssh->term_width, ssh->term_height);
    if (ssh->eof_needed)
	ssh_special(ssh, TS_EOF);

    /*
     * All the initial channel requests are done, so install the default
     * failure handler.
     */
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_SUCCESS] = ssh2_msg_channel_success;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_FAILURE] = ssh2_msg_channel_failure;

    /*
     * Transfer data!
     */
    if (ssh->ldisc)
	ldisc_send(ssh->ldisc, NULL, 0, 0);/* cause ldisc to notice changes */
    if (ssh->mainchan)
	ssh->send_ok = 1;
    while (1) {
	crReturnV;
	s->try_send = FALSE;
	if (pktin) {

	    /*
	     * _All_ the connection-layer packets we expect to
	     * receive are now handled by the dispatch table.
	     * Anything that reaches here must be bogus.
	     */

	    bombout(("Strange packet received: type %d", pktin->type));
	    crStopV;
	} else if (ssh->mainchan) {
	    /*
	     * We have spare data. Add it to the channel buffer.
	     */
	    ssh2_add_channel_data(ssh->mainchan, (char *)in, inlen);
	    s->try_send = TRUE;
	}
	if (s->try_send) {
	    int i;
	    struct ssh_channel *c;
	    /*
	     * Try to send data on all channels if we can.
	     */
	    for (i = 0; NULL != (c = index234(ssh->channels, i)); i++)
		ssh2_try_send_and_unthrottle(ssh, c);
	}
    }

    crFinishV;
}

/*
 * Handlers for SSH-2 messages that might arrive at any moment.
 */
static void ssh2_msg_disconnect(Ssh ssh, struct Packet *pktin)
{
    /* log reason code in disconnect message */
    char *buf, *msg;
    int reason, msglen;

    reason = ssh_pkt_getuint32(pktin);
    ssh_pkt_getstring(pktin, &msg, &msglen);

    if (reason > 0 && reason < lenof(ssh2_disconnect_reasons)) {
	buf = dupprintf("Received disconnect message (%s)",
			ssh2_disconnect_reasons[reason]);
    } else {
	buf = dupprintf("Received disconnect message (unknown"
			" type %d)", reason);
    }
    logevent(buf);
    sfree(buf);
    buf = dupprintf("Disconnection message text: %.*s",
		    msglen, msg);
    logevent(buf);
    bombout(("Server sent disconnect message\ntype %d (%s):\n\"%.*s\"",
	     reason,
	     (reason > 0 && reason < lenof(ssh2_disconnect_reasons)) ?
	     ssh2_disconnect_reasons[reason] : "unknown",
	     msglen, msg));
    sfree(buf);
}

static void ssh2_msg_debug(Ssh ssh, struct Packet *pktin)
{
    /* log the debug message */
    char *msg;
    int msglen;

    /* XXX maybe we should actually take notice of the return value */
    ssh2_pkt_getbool(pktin);
    ssh_pkt_getstring(pktin, &msg, &msglen);

    logeventf(ssh, "Remote debug message: %.*s", msglen, msg);
}

static void ssh2_msg_something_unimplemented(Ssh ssh, struct Packet *pktin)
{
    struct Packet *pktout;
    pktout = ssh2_pkt_init(SSH2_MSG_UNIMPLEMENTED);
    ssh2_pkt_adduint32(pktout, pktin->sequence);
    /*
     * UNIMPLEMENTED messages MUST appear in the same order as the
     * messages they respond to. Hence, never queue them.
     */
    ssh2_pkt_send_noqueue(ssh, pktout);
}

/*
 * Handle the top-level SSH-2 protocol.
 */
static void ssh2_protocol_setup(Ssh ssh)
{
    int i;

    /*
     * Most messages cause SSH2_MSG_UNIMPLEMENTED.
     */
    for (i = 0; i < 256; i++)
	ssh->packet_dispatch[i] = ssh2_msg_something_unimplemented;

    /*
     * Any message we actually understand, we set to NULL so that
     * the coroutines will get it.
     */
    ssh->packet_dispatch[SSH2_MSG_UNIMPLEMENTED] = NULL;
    ssh->packet_dispatch[SSH2_MSG_SERVICE_REQUEST] = NULL;
    ssh->packet_dispatch[SSH2_MSG_SERVICE_ACCEPT] = NULL;
    ssh->packet_dispatch[SSH2_MSG_KEXINIT] = NULL;
    ssh->packet_dispatch[SSH2_MSG_NEWKEYS] = NULL;
    ssh->packet_dispatch[SSH2_MSG_KEXDH_INIT] = NULL;
    ssh->packet_dispatch[SSH2_MSG_KEXDH_REPLY] = NULL;
    /* ssh->packet_dispatch[SSH2_MSG_KEX_DH_GEX_REQUEST] = NULL; duplicate case value */
    /* ssh->packet_dispatch[SSH2_MSG_KEX_DH_GEX_GROUP] = NULL; duplicate case value */
    ssh->packet_dispatch[SSH2_MSG_KEX_DH_GEX_INIT] = NULL;
    ssh->packet_dispatch[SSH2_MSG_KEX_DH_GEX_REPLY] = NULL;
    ssh->packet_dispatch[SSH2_MSG_USERAUTH_REQUEST] = NULL;
    ssh->packet_dispatch[SSH2_MSG_USERAUTH_FAILURE] = NULL;
    ssh->packet_dispatch[SSH2_MSG_USERAUTH_SUCCESS] = NULL;
    ssh->packet_dispatch[SSH2_MSG_USERAUTH_BANNER] = NULL;
    ssh->packet_dispatch[SSH2_MSG_USERAUTH_PK_OK] = NULL;
    /* ssh->packet_dispatch[SSH2_MSG_USERAUTH_PASSWD_CHANGEREQ] = NULL; duplicate case value */
    /* ssh->packet_dispatch[SSH2_MSG_USERAUTH_INFO_REQUEST] = NULL; duplicate case value */
    ssh->packet_dispatch[SSH2_MSG_USERAUTH_INFO_RESPONSE] = NULL;
    ssh->packet_dispatch[SSH2_MSG_GLOBAL_REQUEST] = NULL;
    ssh->packet_dispatch[SSH2_MSG_REQUEST_SUCCESS] = NULL;
    ssh->packet_dispatch[SSH2_MSG_REQUEST_FAILURE] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_OPEN] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_OPEN_CONFIRMATION] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_OPEN_FAILURE] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_WINDOW_ADJUST] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_DATA] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_EXTENDED_DATA] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_EOF] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_CLOSE] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_REQUEST] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_SUCCESS] = NULL;
    ssh->packet_dispatch[SSH2_MSG_CHANNEL_FAILURE] = NULL;

    /*
     * These special message types we install handlers for.
     */
    ssh->packet_dispatch[SSH2_MSG_DISCONNECT] = ssh2_msg_disconnect;
    ssh->packet_dispatch[SSH2_MSG_IGNORE] = ssh_msg_ignore; /* shared with SSH-1 */
    ssh->packet_dispatch[SSH2_MSG_DEBUG] = ssh2_msg_debug;
}

static void ssh2_timer(void *ctx, long now)
{
    Ssh ssh = (Ssh)ctx;

    if (ssh->state == SSH_STATE_CLOSED)
	return;

    if (!ssh->kex_in_progress && ssh->cfg.ssh_rekey_time != 0 &&
	now - ssh->next_rekey >= 0) {
	do_ssh2_transport(ssh, "timeout", -1, NULL);
    }
}

static void ssh2_protocol(Ssh ssh, void *vin, int inlen,
			  struct Packet *pktin)
{
    unsigned char *in = (unsigned char *)vin;
    if (ssh->state == SSH_STATE_CLOSED)
	return;

    if (pktin) {
	ssh->incoming_data_size += pktin->encrypted_len;
	if (!ssh->kex_in_progress &&
	    ssh->max_data_size != 0 &&
	    ssh->incoming_data_size > ssh->max_data_size)
	    do_ssh2_transport(ssh, "too much data received", -1, NULL);
    }

    if (pktin && ssh->packet_dispatch[pktin->type]) {
	ssh->packet_dispatch[pktin->type](ssh, pktin);
	return;
    }

    if (!ssh->protocol_initial_phase_done ||
	(pktin && pktin->type >= 20 && pktin->type < 50)) {
	if (do_ssh2_transport(ssh, in, inlen, pktin) &&
	    !ssh->protocol_initial_phase_done) {
	    ssh->protocol_initial_phase_done = TRUE;
	    /*
	     * Allow authconn to initialise itself.
	     */
	    do_ssh2_authconn(ssh, NULL, 0, NULL);
	}
    } else {
	do_ssh2_authconn(ssh, in, inlen, pktin);
    }
}

/*
 * Called to set up the connection.
 *
 * Returns an error message, or NULL on success.
 */
static const char *ssh_init(void *frontend_handle, void **backend_handle,
			    Config *cfg,
			    char *host, int port, char **realhost, int nodelay,
			    int keepalive)
{
    const char *p;
    Ssh ssh;

    ssh = snew(struct ssh_tag);
    ssh->cfg = *cfg;		       /* STRUCTURE COPY */
    ssh->version = 0;		       /* when not ready yet */
    ssh->s = NULL;
    ssh->cipher = NULL;
    ssh->v1_cipher_ctx = NULL;
    ssh->crcda_ctx = NULL;
    ssh->cscipher = NULL;
    ssh->cs_cipher_ctx = NULL;
    ssh->sccipher = NULL;
    ssh->sc_cipher_ctx = NULL;
    ssh->csmac = NULL;
    ssh->cs_mac_ctx = NULL;
    ssh->scmac = NULL;
    ssh->sc_mac_ctx = NULL;
    ssh->cscomp = NULL;
    ssh->cs_comp_ctx = NULL;
    ssh->sccomp = NULL;
    ssh->sc_comp_ctx = NULL;
    ssh->kex = NULL;
    ssh->kex_ctx = NULL;
    ssh->hostkey = NULL;
    ssh->exitcode = -1;
    ssh->close_expected = FALSE;
    ssh->clean_exit = FALSE;
    ssh->state = SSH_STATE_PREPACKET;
    ssh->size_needed = FALSE;
    ssh->eof_needed = FALSE;
    ssh->ldisc = NULL;
    ssh->logctx = NULL;
    ssh->deferred_send_data = NULL;
    ssh->deferred_len = 0;
    ssh->deferred_size = 0;
    ssh->fallback_cmd = 0;
    ssh->pkt_kctx = SSH2_PKTCTX_NOKEX;
    ssh->pkt_actx = SSH2_PKTCTX_NOAUTH;
    ssh->x11disp = NULL;
    ssh->v1_compressing = FALSE;
    ssh->v2_outgoing_sequence = 0;
    ssh->ssh1_rdpkt_crstate = 0;
    ssh->ssh2_rdpkt_crstate = 0;
    ssh->do_ssh_init_crstate = 0;
    ssh->ssh_gotdata_crstate = 0;
    ssh->do_ssh1_connection_crstate = 0;
    ssh->do_ssh1_login_crstate = 0;
    ssh->do_ssh2_transport_crstate = 0;
    ssh->do_ssh2_authconn_crstate = 0;
    ssh->do_ssh_init_state = NULL;
    ssh->do_ssh1_login_state = NULL;
    ssh->do_ssh2_transport_state = NULL;
    ssh->do_ssh2_authconn_state = NULL;
    ssh->v_c = NULL;
    ssh->v_s = NULL;
    ssh->mainchan = NULL;
    ssh->throttled_all = 0;
    ssh->v1_stdout_throttling = 0;
    ssh->queue = NULL;
    ssh->queuelen = ssh->queuesize = 0;
    ssh->queueing = FALSE;
    ssh->qhead = ssh->qtail = NULL;
    ssh->deferred_rekey_reason = NULL;
    bufchain_init(&ssh->queued_incoming_data);
    ssh->frozen = FALSE;

    *backend_handle = ssh;

#ifdef MSCRYPTOAPI
    if (crypto_startup() == 0)
	return "Microsoft high encryption pack not installed!";
#endif

    ssh->frontend = frontend_handle;
    ssh->term_width = ssh->cfg.width;
    ssh->term_height = ssh->cfg.height;

    ssh->channels = NULL;
    ssh->rportfwds = NULL;
    ssh->portfwds = NULL;

    ssh->send_ok = 0;
    ssh->editing = 0;
    ssh->echoing = 0;
    ssh->conn_throttle_count = 0;
    ssh->overall_bufsize = 0;
    ssh->fallback_cmd = 0;

    ssh->protocol = NULL;

    ssh->protocol_initial_phase_done = FALSE;

    ssh->pinger = NULL;

    ssh->incoming_data_size = ssh->outgoing_data_size =
	ssh->deferred_data_size = 0L;
    ssh->max_data_size = parse_blocksize(ssh->cfg.ssh_rekey_data);
    ssh->kex_in_progress = FALSE;

#ifndef NO_GSSAPI
    ssh->gsslibs = NULL;
#endif

    p = connect_to_host(ssh, host, port, realhost, nodelay, keepalive);
    if (p != NULL)
	return p;

    random_ref();

    return NULL;
}

static void ssh_free(void *handle)
{
    Ssh ssh = (Ssh) handle;
    struct ssh_channel *c;
    struct ssh_rportfwd *pf;

    if (ssh->v1_cipher_ctx)
	ssh->cipher->free_context(ssh->v1_cipher_ctx);
    if (ssh->cs_cipher_ctx)
	ssh->cscipher->free_context(ssh->cs_cipher_ctx);
    if (ssh->sc_cipher_ctx)
	ssh->sccipher->free_context(ssh->sc_cipher_ctx);
    if (ssh->cs_mac_ctx)
	ssh->csmac->free_context(ssh->cs_mac_ctx);
    if (ssh->sc_mac_ctx)
	ssh->scmac->free_context(ssh->sc_mac_ctx);
    if (ssh->cs_comp_ctx) {
	if (ssh->cscomp)
	    ssh->cscomp->compress_cleanup(ssh->cs_comp_ctx);
	else
	    zlib_compress_cleanup(ssh->cs_comp_ctx);
    }
    if (ssh->sc_comp_ctx) {
	if (ssh->sccomp)
	    ssh->sccomp->decompress_cleanup(ssh->sc_comp_ctx);
	else
	    zlib_decompress_cleanup(ssh->sc_comp_ctx);
    }
    if (ssh->kex_ctx)
	dh_cleanup(ssh->kex_ctx);
    sfree(ssh->savedhost);

    while (ssh->queuelen-- > 0)
	ssh_free_packet(ssh->queue[ssh->queuelen]);
    sfree(ssh->queue);

    while (ssh->qhead) {
	struct queued_handler *qh = ssh->qhead;
	ssh->qhead = qh->next;
	sfree(ssh->qhead);
    }
    ssh->qhead = ssh->qtail = NULL;

    if (ssh->channels) {
	while ((c = delpos234(ssh->channels, 0)) != NULL) {
	    switch (c->type) {
	      case CHAN_X11:
		if (c->u.x11.s != NULL)
		    x11_close(c->u.x11.s);
		break;
	      case CHAN_SOCKDATA:
	      case CHAN_SOCKDATA_DORMANT:
		if (c->u.pfd.s != NULL)
		    pfd_close(c->u.pfd.s);
		break;
	    }
	    sfree(c);
	}
	freetree234(ssh->channels);
	ssh->channels = NULL;
    }

    if (ssh->rportfwds) {
	while ((pf = delpos234(ssh->rportfwds, 0)) != NULL)
	    free_rportfwd(pf);
	freetree234(ssh->rportfwds);
	ssh->rportfwds = NULL;
    }
    sfree(ssh->deferred_send_data);
    if (ssh->x11disp)
	x11_free_display(ssh->x11disp);
    sfree(ssh->do_ssh_init_state);
    sfree(ssh->do_ssh1_login_state);
    sfree(ssh->do_ssh2_transport_state);
    sfree(ssh->do_ssh2_authconn_state);
    sfree(ssh->v_c);
    sfree(ssh->v_s);
    sfree(ssh->fullhostname);
    if (ssh->crcda_ctx) {
	crcda_free_context(ssh->crcda_ctx);
	ssh->crcda_ctx = NULL;
    }
    if (ssh->s)
	ssh_do_close(ssh, TRUE);
    expire_timer_context(ssh);
    if (ssh->pinger)
	pinger_free(ssh->pinger);
    bufchain_clear(&ssh->queued_incoming_data);
#ifndef NO_GSSAPI
    if (ssh->gsslibs)
	ssh_gss_cleanup(ssh->gsslibs);
#endif
    sfree(ssh);

    random_unref();
}

/*
 * Reconfigure the SSH backend.
 */
static void ssh_reconfig(void *handle, Config *cfg)
{
    Ssh ssh = (Ssh) handle;
    char *rekeying = NULL, rekey_mandatory = FALSE;
    unsigned long old_max_data_size;

    pinger_reconfig(ssh->pinger, &ssh->cfg, cfg);
    if (ssh->portfwds)
	ssh_setup_portfwd(ssh, cfg);

    if (ssh->cfg.ssh_rekey_time != cfg->ssh_rekey_time &&
	cfg->ssh_rekey_time != 0) {
	long new_next = ssh->last_rekey + cfg->ssh_rekey_time*60*TICKSPERSEC;
	long now = GETTICKCOUNT();

	if (new_next - now < 0) {
	    rekeying = "timeout shortened";
	} else {
	    ssh->next_rekey = schedule_timer(new_next - now, ssh2_timer, ssh);
	}
    }

    old_max_data_size = ssh->max_data_size;
    ssh->max_data_size = parse_blocksize(cfg->ssh_rekey_data);
    if (old_max_data_size != ssh->max_data_size &&
	ssh->max_data_size != 0) {
	if (ssh->outgoing_data_size > ssh->max_data_size ||
	    ssh->incoming_data_size > ssh->max_data_size)
	    rekeying = "data limit lowered";
    }

    if (ssh->cfg.compression != cfg->compression) {
	rekeying = "compression setting changed";
	rekey_mandatory = TRUE;
    }

    if (ssh->cfg.ssh2_des_cbc != cfg->ssh2_des_cbc ||
	memcmp(ssh->cfg.ssh_cipherlist, cfg->ssh_cipherlist,
	       sizeof(ssh->cfg.ssh_cipherlist))) {
	rekeying = "cipher settings changed";
	rekey_mandatory = TRUE;
    }

    ssh->cfg = *cfg;		       /* STRUCTURE COPY */

    if (rekeying) {
	if (!ssh->kex_in_progress) {
	    do_ssh2_transport(ssh, rekeying, -1, NULL);
	} else if (rekey_mandatory) {
	    ssh->deferred_rekey_reason = rekeying;
	}
    }
}

/*
 * Called to send data down the SSH connection.
 */
static int ssh_send(void *handle, char *buf, int len)
{
    Ssh ssh = (Ssh) handle;

    if (ssh == NULL || ssh->s == NULL || ssh->protocol == NULL)
	return 0;

    ssh->protocol(ssh, (unsigned char *)buf, len, 0);

    return ssh_sendbuffer(ssh);
}

/*
 * Called to query the current amount of buffered stdin data.
 */
static int ssh_sendbuffer(void *handle)
{
    Ssh ssh = (Ssh) handle;
    int override_value;

    if (ssh == NULL || ssh->s == NULL || ssh->protocol == NULL)
	return 0;

    /*
     * If the SSH socket itself has backed up, add the total backup
     * size on that to any individual buffer on the stdin channel.
     */
    override_value = 0;
    if (ssh->throttled_all)
	override_value = ssh->overall_bufsize;

    if (ssh->version == 1) {
	return override_value;
    } else if (ssh->version == 2) {
	if (!ssh->mainchan || ssh->mainchan->closes > 0)
	    return override_value;
	else
	    return (override_value +
		    bufchain_size(&ssh->mainchan->v.v2.outbuffer));
    }

    return 0;
}

/*
 * Called to set the size of the window from SSH's POV.
 */
static void ssh_size(void *handle, int width, int height)
{
    Ssh ssh = (Ssh) handle;
    struct Packet *pktout;

    ssh->term_width = width;
    ssh->term_height = height;

    switch (ssh->state) {
      case SSH_STATE_BEFORE_SIZE:
      case SSH_STATE_PREPACKET:
      case SSH_STATE_CLOSED:
	break;			       /* do nothing */
      case SSH_STATE_INTERMED:
	ssh->size_needed = TRUE;       /* buffer for later */
	break;
      case SSH_STATE_SESSION:
	if (!ssh->cfg.nopty) {
	    if (ssh->version == 1) {
		send_packet(ssh, SSH1_CMSG_WINDOW_SIZE,
			    PKT_INT, ssh->term_height,
			    PKT_INT, ssh->term_width,
			    PKT_INT, 0, PKT_INT, 0, PKT_END);
	    } else if (ssh->mainchan) {
		pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
		ssh2_pkt_adduint32(pktout, ssh->mainchan->remoteid);
		ssh2_pkt_addstring(pktout, "window-change");
		ssh2_pkt_addbool(pktout, 0);
		ssh2_pkt_adduint32(pktout, ssh->term_width);
		ssh2_pkt_adduint32(pktout, ssh->term_height);
		ssh2_pkt_adduint32(pktout, 0);
		ssh2_pkt_adduint32(pktout, 0);
		ssh2_pkt_send(ssh, pktout);
	    }
	}
	break;
    }
}

/*
 * Return a list of the special codes that make sense in this
 * protocol.
 */
static const struct telnet_special *ssh_get_specials(void *handle)
{
    static const struct telnet_special ssh1_ignore_special[] = {
	{"IGNORE message", TS_NOP}
    };
    static const struct telnet_special ssh2_ignore_special[] = {
	{"IGNORE message", TS_NOP},
    };
    static const struct telnet_special ssh2_rekey_special[] = {
	{"Repeat key exchange", TS_REKEY},
    };
    static const struct telnet_special ssh2_session_specials[] = {
	{NULL, TS_SEP},
	{"Break", TS_BRK},
	/* These are the signal names defined by RFC 4254.
	 * They include all the ISO C signals, but are a subset of the POSIX
	 * required signals. */
	{"SIGINT (Interrupt)", TS_SIGINT},
	{"SIGTERM (Terminate)", TS_SIGTERM},
	{"SIGKILL (Kill)", TS_SIGKILL},
	{"SIGQUIT (Quit)", TS_SIGQUIT},
	{"SIGHUP (Hangup)", TS_SIGHUP},
	{"More signals", TS_SUBMENU},
	  {"SIGABRT", TS_SIGABRT}, {"SIGALRM", TS_SIGALRM},
	  {"SIGFPE",  TS_SIGFPE},  {"SIGILL",  TS_SIGILL},
	  {"SIGPIPE", TS_SIGPIPE}, {"SIGSEGV", TS_SIGSEGV},
	  {"SIGUSR1", TS_SIGUSR1}, {"SIGUSR2", TS_SIGUSR2},
	{NULL, TS_EXITMENU}
    };
    static const struct telnet_special specials_end[] = {
	{NULL, TS_EXITMENU}
    };
    /* XXX review this length for any changes: */
    static struct telnet_special ssh_specials[lenof(ssh2_ignore_special) +
					      lenof(ssh2_rekey_special) +
					      lenof(ssh2_session_specials) +
					      lenof(specials_end)];
    Ssh ssh = (Ssh) handle;
    int i = 0;
#define ADD_SPECIALS(name) \
    do { \
	assert((i + lenof(name)) <= lenof(ssh_specials)); \
	memcpy(&ssh_specials[i], name, sizeof name); \
	i += lenof(name); \
    } while(0)

    if (ssh->version == 1) {
	/* Don't bother offering IGNORE if we've decided the remote
	 * won't cope with it, since we wouldn't bother sending it if
	 * asked anyway. */
	if (!(ssh->remote_bugs & BUG_CHOKES_ON_SSH1_IGNORE))
	    ADD_SPECIALS(ssh1_ignore_special);
    } else if (ssh->version == 2) {
	if (!(ssh->remote_bugs & BUG_CHOKES_ON_SSH2_IGNORE))
	    ADD_SPECIALS(ssh2_ignore_special);
	if (!(ssh->remote_bugs & BUG_SSH2_REKEY))
	    ADD_SPECIALS(ssh2_rekey_special);
	if (ssh->mainchan)
	    ADD_SPECIALS(ssh2_session_specials);
    } /* else we're not ready yet */

    if (i) {
	ADD_SPECIALS(specials_end);
	return ssh_specials;
    } else {
	return NULL;
    }
#undef ADD_SPECIALS
}

/*
 * Send special codes. TS_EOF is useful for `plink', so you
 * can send an EOF and collect resulting output (e.g. `plink
 * hostname sort').
 */
static void ssh_special(void *handle, Telnet_Special code)
{
    Ssh ssh = (Ssh) handle;
    struct Packet *pktout;

    if (code == TS_EOF) {
	if (ssh->state != SSH_STATE_SESSION) {
	    /*
	     * Buffer the EOF in case we are pre-SESSION, so we can
	     * send it as soon as we reach SESSION.
	     */
	    if (code == TS_EOF)
		ssh->eof_needed = TRUE;
	    return;
	}
	if (ssh->version == 1) {
	    send_packet(ssh, SSH1_CMSG_EOF, PKT_END);
	} else if (ssh->mainchan) {
	    struct Packet *pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_EOF);
	    ssh2_pkt_adduint32(pktout, ssh->mainchan->remoteid);
	    ssh2_pkt_send(ssh, pktout);
            ssh->send_ok = 0;          /* now stop trying to read from stdin */
	}
	logevent("Sent EOF message");
    } else if (code == TS_PING || code == TS_NOP) {
	if (ssh->state == SSH_STATE_CLOSED
	    || ssh->state == SSH_STATE_PREPACKET) return;
	if (ssh->version == 1) {
	    if (!(ssh->remote_bugs & BUG_CHOKES_ON_SSH1_IGNORE))
		send_packet(ssh, SSH1_MSG_IGNORE, PKT_STR, "", PKT_END);
	} else {
	    if (!(ssh->remote_bugs & BUG_CHOKES_ON_SSH2_IGNORE)) {
		pktout = ssh2_pkt_init(SSH2_MSG_IGNORE);
		ssh2_pkt_addstring_start(pktout);
		ssh2_pkt_send_noqueue(ssh, pktout);
	    }
	}
    } else if (code == TS_REKEY) {
	if (!ssh->kex_in_progress && ssh->version == 2) {
	    do_ssh2_transport(ssh, "at user request", -1, NULL);
	}
    } else if (code == TS_BRK) {
	if (ssh->state == SSH_STATE_CLOSED
	    || ssh->state == SSH_STATE_PREPACKET) return;
	if (ssh->version == 1) {
	    logevent("Unable to send BREAK signal in SSH-1");
	} else if (ssh->mainchan) {
	    pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
	    ssh2_pkt_adduint32(pktout, ssh->mainchan->remoteid);
	    ssh2_pkt_addstring(pktout, "break");
	    ssh2_pkt_addbool(pktout, 0);
	    ssh2_pkt_adduint32(pktout, 0);   /* default break length */
	    ssh2_pkt_send(ssh, pktout);
	}
    } else {
	/* Is is a POSIX signal? */
	char *signame = NULL;
	if (code == TS_SIGABRT) signame = "ABRT";
	if (code == TS_SIGALRM) signame = "ALRM";
	if (code == TS_SIGFPE)  signame = "FPE";
	if (code == TS_SIGHUP)  signame = "HUP";
	if (code == TS_SIGILL)  signame = "ILL";
	if (code == TS_SIGINT)  signame = "INT";
	if (code == TS_SIGKILL) signame = "KILL";
	if (code == TS_SIGPIPE) signame = "PIPE";
	if (code == TS_SIGQUIT) signame = "QUIT";
	if (code == TS_SIGSEGV) signame = "SEGV";
	if (code == TS_SIGTERM) signame = "TERM";
	if (code == TS_SIGUSR1) signame = "USR1";
	if (code == TS_SIGUSR2) signame = "USR2";
	/* The SSH-2 protocol does in principle support arbitrary named
	 * signals, including signame@domain, but we don't support those. */
	if (signame) {
	    /* It's a signal. */
	    if (ssh->version == 2 && ssh->mainchan) {
		pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_REQUEST);
		ssh2_pkt_adduint32(pktout, ssh->mainchan->remoteid);
		ssh2_pkt_addstring(pktout, "signal");
		ssh2_pkt_addbool(pktout, 0);
		ssh2_pkt_addstring(pktout, signame);
		ssh2_pkt_send(ssh, pktout);
		logeventf(ssh, "Sent signal SIG%s", signame);
	    }
	} else {
	    /* Never heard of it. Do nothing */
	}
    }
}

void *new_sock_channel(void *handle, Socket s)
{
    Ssh ssh = (Ssh) handle;
    struct ssh_channel *c;
    c = snew(struct ssh_channel);

    c->ssh = ssh;
    ssh2_channel_init(c);
    c->halfopen = TRUE;
    c->type = CHAN_SOCKDATA_DORMANT;/* identify channel type */
    c->u.pfd.s = s;
    add234(ssh->channels, c);
    return c;
}

/*
 * This is called when stdout/stderr (the entity to which
 * from_backend sends data) manages to clear some backlog.
 */
static void ssh_unthrottle(void *handle, int bufsize)
{
    Ssh ssh = (Ssh) handle;
    int buflimit;

    if (ssh->version == 1) {
	if (ssh->v1_stdout_throttling && bufsize < SSH1_BUFFER_LIMIT) {
	    ssh->v1_stdout_throttling = 0;
	    ssh_throttle_conn(ssh, -1);
	}
    } else {
	if (ssh->mainchan) {
	    ssh2_set_window(ssh->mainchan,
			    bufsize < ssh->mainchan->v.v2.locmaxwin ?
			    ssh->mainchan->v.v2.locmaxwin - bufsize : 0);
	    if (ssh->cfg.ssh_simple)
		buflimit = 0;
	    else
		buflimit = ssh->mainchan->v.v2.locmaxwin;
	    if (ssh->mainchan->throttling_conn && bufsize <= buflimit) {
		ssh->mainchan->throttling_conn = 0;
		ssh_throttle_conn(ssh, -1);
	    }
	}
    }
}

void ssh_send_port_open(void *channel, char *hostname, int port, char *org)
{
    struct ssh_channel *c = (struct ssh_channel *)channel;
    Ssh ssh = c->ssh;
    struct Packet *pktout;

    logeventf(ssh, "Opening forwarded connection to %s:%d", hostname, port);

    if (ssh->version == 1) {
	send_packet(ssh, SSH1_MSG_PORT_OPEN,
		    PKT_INT, c->localid,
		    PKT_STR, hostname,
		    PKT_INT, port,
		    /* PKT_STR, <org:orgport>, */
		    PKT_END);
    } else {
	pktout = ssh2_pkt_init(SSH2_MSG_CHANNEL_OPEN);
	ssh2_pkt_addstring(pktout, "direct-tcpip");
	ssh2_pkt_adduint32(pktout, c->localid);
	ssh2_pkt_adduint32(pktout, c->v.v2.locwindow);/* our window size */
	ssh2_pkt_adduint32(pktout, OUR_V2_MAXPKT);      /* our max pkt size */
	ssh2_pkt_addstring(pktout, hostname);
	ssh2_pkt_adduint32(pktout, port);
	/*
	 * We make up values for the originator data; partly it's
	 * too much hassle to keep track, and partly I'm not
	 * convinced the server should be told details like that
	 * about my local network configuration.
	 * The "originator IP address" is syntactically a numeric
	 * IP address, and some servers (e.g., Tectia) get upset
	 * if it doesn't match this syntax.
	 */
	ssh2_pkt_addstring(pktout, "0.0.0.0");
	ssh2_pkt_adduint32(pktout, 0);
	ssh2_pkt_send(ssh, pktout);
    }
}

static int ssh_connected(void *handle)
{
    Ssh ssh = (Ssh) handle;
    return ssh->s != NULL;
}

static int ssh_sendok(void *handle)
{
    Ssh ssh = (Ssh) handle;
    return ssh->send_ok;
}

static int ssh_ldisc(void *handle, int option)
{
    Ssh ssh = (Ssh) handle;
    if (option == LD_ECHO)
	return ssh->echoing;
    if (option == LD_EDIT)
	return ssh->editing;
    return FALSE;
}

static void ssh_provide_ldisc(void *handle, void *ldisc)
{
    Ssh ssh = (Ssh) handle;
    ssh->ldisc = ldisc;
}

static void ssh_provide_logctx(void *handle, void *logctx)
{
    Ssh ssh = (Ssh) handle;
    ssh->logctx = logctx;
}

static int ssh_return_exitcode(void *handle)
{
    Ssh ssh = (Ssh) handle;
    if (ssh->s != NULL)
        return -1;
    else
        return (ssh->exitcode >= 0 ? ssh->exitcode : INT_MAX);
}

/*
 * cfg_info for SSH is the currently running version of the
 * protocol. (1 for 1; 2 for 2; 0 for not-decided-yet.)
 */
static int ssh_cfg_info(void *handle)
{
    Ssh ssh = (Ssh) handle;
    return ssh->version;
}

/*
 * Gross hack: pscp will try to start SFTP but fall back to scp1 if
 * that fails. This variable is the means by which scp.c can reach
 * into the SSH code and find out which one it got.
 */
extern int ssh_fallback_cmd(void *handle)
{
    Ssh ssh = (Ssh) handle;
    return ssh->fallback_cmd;
}

Backend ssh_backend = {
    ssh_init,
    ssh_free,
    ssh_reconfig,
    ssh_send,
    ssh_sendbuffer,
    ssh_size,
    ssh_special,
    ssh_get_specials,
    ssh_connected,
    ssh_return_exitcode,
    ssh_sendok,
    ssh_ldisc,
    ssh_provide_ldisc,
    ssh_provide_logctx,
    ssh_unthrottle,
    ssh_cfg_info,
    "ssh",
    PROT_SSH,
    22
};
