import { Card, CardContent, CardHeader } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import { CalendarDays, Clock, User, ExternalLink, AlertTriangle, MapPin, Users, Video } from "lucide-react"
import { format } from "date-fns"
import { zhCN } from "date-fns/locale"
import Link from "next/link"
import Image from "next/image"
import type { Announcement, Post, LearnResource, Meeting } from "@/lib/contentlayer-mock"

interface AnnouncementCardProps {
  announcement: Announcement
}

export function AnnouncementCard({ announcement }: AnnouncementCardProps) {
  return (
    <Card className={`transition-all hover:shadow-md ${announcement.isExpired ? "opacity-75" : ""}`}>
      <CardHeader className="pb-3">
        <div className="flex items-start justify-between gap-4">
          <div className="flex-1">
            <div className="flex flex-wrap gap-2 mb-3">
              {announcement.tags.map((tag) => (
                <Badge key={tag} variant="secondary" className="text-xs">
                  {tag}
                </Badge>
              ))}
              {announcement.priority > 0 && (
                <Badge variant="destructive" className="text-xs">
                  置顶
                </Badge>
              )}
              {announcement.isExpired && (
                <Badge variant="outline" className="text-xs text-muted-foreground">
                  <AlertTriangle className="h-3 w-3 mr-1" />
                  已过期
                </Badge>
              )}
            </div>
            <h3 className="text-lg font-semibold text-balance leading-tight mb-2">
              <Link href={announcement.url} className="hover:text-primary transition-colors">
                {announcement.title}
              </Link>
            </h3>
          </div>
        </div>
      </CardHeader>
      <CardContent>
        <p className="text-muted-foreground text-sm text-pretty mb-4 line-clamp-3">{announcement.summary}</p>
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4 text-xs text-muted-foreground">
            <div className="flex items-center gap-1">
              <CalendarDays className="h-3 w-3" />
              <span>{format(new Date(announcement.publishedAt), "MM/dd", { locale: zhCN })}</span>
            </div>
            {announcement.expiresAt && (
              <div className="flex items-center gap-1">
                <Clock className="h-3 w-3" />
                <span>截止 {format(new Date(announcement.expiresAt), "MM/dd", { locale: zhCN })}</span>
              </div>
            )}
          </div>
          <Button asChild variant="ghost" size="sm">
            <Link href={announcement.url}>查看详情</Link>
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}

interface PostCardProps {
  post: Post
}

export function PostCard({ post }: PostCardProps) {
  return (
    <Card className="transition-all hover:shadow-md">
      {post.coverImage && (
        <div className="aspect-video overflow-hidden rounded-t-lg">
          <Image
            src={post.coverImage || "/placeholder.svg"}
            alt={post.title}
            width={400}
            height={225}
            className="w-full h-full object-cover transition-transform hover:scale-105"
            unoptimized
          />
        </div>
      )}
      <CardHeader className="pb-3">
        <div className="flex flex-wrap gap-2 mb-3">
          {post.tags.slice(0, 3).map((tag) => (
            <Badge key={tag} variant="secondary" className="text-xs">
              {tag}
            </Badge>
          ))}
          {post.tags.length > 3 && (
            <Badge variant="outline" className="text-xs">
              +{post.tags.length - 3}
            </Badge>
          )}
        </div>
        <h3 className="text-lg font-semibold text-balance leading-tight mb-2">
          <Link href={post.url} className="hover:text-primary transition-colors">
            {post.title}
          </Link>
        </h3>
      </CardHeader>
      <CardContent>
        <p className="text-muted-foreground text-sm text-pretty mb-4 line-clamp-3">{post.excerpt}</p>
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4 text-xs text-muted-foreground">
            <div className="flex items-center gap-1">
              <User className="h-3 w-3" />
              <span>{post.authors[0]}</span>
              {post.authors.length > 1 && <span>等</span>}
            </div>
            <div className="flex items-center gap-1">
              <CalendarDays className="h-3 w-3" />
              <span>{format(new Date(post.publishedAt), "MM/dd", { locale: zhCN })}</span>
            </div>
          </div>
          <Button asChild variant="ghost" size="sm">
            <Link href={post.url}>阅读全文</Link>
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}

interface LearnResourceCardProps {
  resource: LearnResource
}

export function LearnResourceCard({ resource }: LearnResourceCardProps) {
  const levelColors = {
    入门: "bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200",
    进阶: "bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200",
    实战: "bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-200",
  }

  return (
    <Card className="transition-all hover:shadow-md">
      <CardHeader className="pb-3">
        <div className="flex items-start justify-between gap-4">
          <div className="flex-1">
            <div className="flex flex-wrap gap-2 mb-3">
              <Badge variant="outline" className="text-xs">
                {resource.type}
              </Badge>
              <Badge className={`text-xs ${levelColors[resource.level]}`}>{resource.level}</Badge>
              <Badge variant="secondary" className="text-xs">
                {resource.track}
              </Badge>
            </div>
            <h3 className="text-lg font-semibold text-balance leading-tight mb-2">
              <Link href={resource.url} className="hover:text-primary transition-colors">
                {resource.title}
              </Link>
            </h3>
          </div>
          {resource.externalUrl && <ExternalLink className="h-4 w-4 text-muted-foreground flex-shrink-0" />}
        </div>
      </CardHeader>
      <CardContent>
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4 text-xs text-muted-foreground">
            <div className="flex items-center gap-1">
              <CalendarDays className="h-3 w-3" />
              <span>{format(new Date(resource.updatedAt), "MM/dd", { locale: zhCN })}</span>
            </div>
            {resource.duration && (
              <div className="flex items-center gap-1">
                <Clock className="h-3 w-3" />
                <span>{resource.duration}</span>
              </div>
            )}
          </div>
          <Button asChild variant="ghost" size="sm">
            <Link href={resource.url}>{resource.externalUrl ? "访问资源" : "查看详情"}</Link>
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}

interface MeetingCardProps {
  meeting: Meeting
}

export function MeetingCard({ meeting }: MeetingCardProps) {
  const statusColors = {
    即将开始: "bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200",
    进行中: "bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200",
    已结束: "bg-gray-100 text-gray-800 dark:bg-gray-900 dark:text-gray-200",
  }

  const typeIcons = {
    线上: Video,
    线下: MapPin,
    混合: Users,
  }

  const TypeIcon = typeIcons[meeting.type]

  return (
    <Card className={`transition-all hover:shadow-md ${meeting.status === "已结束" ? "opacity-75" : ""}`}>
      <CardHeader className="pb-3">
        <div className="flex items-start justify-between gap-4">
          <div className="flex-1">
            <div className="flex flex-wrap gap-2 mb-3">
              {meeting.tags.map((tag) => (
                <Badge key={tag} variant="secondary" className="text-xs">
                  {tag}
                </Badge>
              ))}
              <Badge className={`text-xs ${statusColors[meeting.status]}`}>{meeting.status}</Badge>
            </div>
            <h3 className="text-lg font-semibold text-balance leading-tight mb-2">
              <Link href={meeting.url} className="hover:text-primary transition-colors">
                {meeting.title}
              </Link>
            </h3>
          </div>
        </div>
      </CardHeader>
      <CardContent>
        <p className="text-muted-foreground text-sm text-pretty mb-4 line-clamp-3">{meeting.summary}</p>

        <div className="space-y-2 mb-4">
          <div className="flex items-center gap-2 text-sm">
            <CalendarDays className="h-4 w-4 text-muted-foreground" />
            <span>{format(new Date(meeting.date), "yyyy年MM月dd日", { locale: zhCN })}</span>
            <Clock className="h-4 w-4 text-muted-foreground ml-2" />
            <span>{meeting.time}</span>
          </div>

          <div className="flex items-center gap-2 text-sm">
            <TypeIcon className="h-4 w-4 text-muted-foreground" />
            <span>{meeting.location}</span>
            <Badge variant="outline" className="text-xs ml-2">
              {meeting.type}
            </Badge>
          </div>

          {meeting.speaker && (
            <div className="flex items-center gap-2 text-sm">
              <User className="h-4 w-4 text-muted-foreground" />
              <span>主讲：{meeting.speaker}</span>
            </div>
          )}

          {meeting.capacity && (
            <div className="flex items-center gap-2 text-sm">
              <Users className="h-4 w-4 text-muted-foreground" />
              <span>
                {meeting.registered || 0}/{meeting.capacity} 人
                {meeting.registered && meeting.capacity && (
                  <span className="text-muted-foreground ml-1">
                    ({Math.round((meeting.registered / meeting.capacity) * 100)}% 已报名)
                  </span>
                )}
              </span>
            </div>
          )}
        </div>

        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            {meeting.meetingUrl && meeting.status !== "已结束" && (
              <Button asChild variant="outline" size="sm">
                <Link href={meeting.meetingUrl} target="_blank" rel="noopener noreferrer">
                  <Video className="h-3 w-3 mr-1" />
                  加入会议
                </Link>
              </Button>
            )}
          </div>
          <Button asChild variant="ghost" size="sm">
            <Link href={meeting.url}>查看详情</Link>
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}
